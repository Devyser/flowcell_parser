import re
import os
import csv
import xml.etree.ElementTree as ET
import logging
import glob
import json
import pandas as pd
from datetime import datetime
from collections import OrderedDict
from bs4 import BeautifulSoup #html parser
from interop.core import imaging
from interop.py_interop_run import xml_file_not_found_exception


class RunParser(object):
    """Parses an Illumina run folder. It generates data for statusdb
    notable attributes :

    :RunInfoParser runinfo: see RunInfo
    :RunParametersParser runparameters: see RunParametersParser
    :SampleSheetParser samplesheet: see SampleSheetParser
    :LaneBarcodeParser lanebarcodes: see LaneBarcodeParser
    """

    def __init__(self, path):
        if os.path.exists(path):
            self.log = logging.getLogger(__name__)
            self.path = path
            self.parse()
            self.create_db_obj()
        else:
            raise os.error(f"Flowcell cannot be found at {path}")

    def parse(self, demultiplexingDir = 'Demultiplexing'):
        """Tries to parse as many files as possible from a run folder"""
        # When demultiplexing with BCL Convert, e.g. used for MiSeq i100, the old reports can be produced with '--output-legacy-stats true',
        # but they will end up in a different location. Also, BCL Convert uses a new sample sheet format.
        pattern = r'(\d{6,8})_([ST-]*\w+\d+)_\d+_([AB]?)([A-Z0-9\-]+)'
        pattern_match = re.match(pattern, os.path.basename(os.path.abspath(self.path)))
        if pattern_match.group(2).startswith('SL'):
            demultiplexingDir = os.path.join(demultiplexingDir, 'Reports', 'legacy')
            samplesheet_parser = SampleSheetV2Parser
        else:
            samplesheet_parser = SampleSheetParser

        # Start with the files existing before demultiplexing, located directly in run folder
        rinfo_path = os.path.join(self.path, 'RunInfo.xml')
        rpar_path = os.path.join(self.path, 'runParameters.xml')
        if not os.path.exists(rpar_path):
            rpar_path = os.path.join(self.path, 'RunParameters.xml')
        ss_path = os.path.join(self.path, 'SampleSheet.csv')
        cycle_times_log = os.path.join(self.path, 'Logs', 'CycleTimes.txt')

        try:
            self.runinfo = RunInfoParser(rinfo_path)
        except OSError as e:
            self.log.info(str(e))
            self.runinfo = None
        try:
            self.runparameters = RunParametersParser(rpar_path)
        except OSError as e:
            self.log.info(str(e))
            self.runparameters = None
        try:
            self.samplesheet = samplesheet_parser(ss_path)
        except OSError as e:
            self.log.info(str(e))
            self.samplesheet = None
        try:
            self.interop_data = InterOpParser(self.path)
        except (OSError, xml_file_not_found_exception):
            self.log.info(str(e))
            self.interop_data = None

        # Continue with files generate post-demultiplexing and could thus potentially be replaced by reading from stats.json
        try:
            if self.runinfo:
                fc_name = self.runinfo.data.get('Flowcell','')
            elif len(pattern_match.group(1)) == 8:
                # Try to continue for an iSeq flow cell
                fc_name = pattern_match.group(3)
            else:
                # Try to continue for a MiSeq/MiniSeq flow cell
                fc_name = pattern_match.group(4)
        except IndexError:
            self.log.error("Was not able to identify flow cell name. Flow cell directory has an unexpected pattern")
            fc_name = ''

        lb_path = os.path.join(self.path, demultiplexingDir, 'Reports', 'html', fc_name, 'all', 'all', 'all', 'laneBarcode.html')
        ln_path = os.path.join(self.path, demultiplexingDir, 'Reports', 'html', fc_name, 'all', 'all', 'all', 'lane.html')
        undeterminedStatsFolder = os.path.join(self.path, demultiplexingDir, 'Stats')
        json_path = os.path.join(self.path, demultiplexingDir, 'Stats', 'Stats.json')

        try:
            self.lanebarcodes = LaneBarcodeParser(lb_path)
        except OSError as e:
            self.log.info(str(e))
            self.lanebarcodes = None
        try:
            self.lanes = LaneBarcodeParser(ln_path)
        except OSError as e:
            self.log.info(str(e))
            self.lanes = None
        try:
            self.undet = DemuxSummaryParser(undeterminedStatsFolder)
        except OSError as e:
            self.log.info(str(e))
            self.undet = None
        try:
            self.time_cycles = CycleTimesParser(cycle_times_log)
        except OSError as e:
            self.log.info(str(e))
            self.time_cycles = None
        try:
            self.json_stats = StatsParser(json_path)
        except OSError as e:
            self.log.info(str(e))
            self.json_stats = None

    def create_db_obj(self):        
        self.obj = {}
   
        bits = os.path.basename(os.path.abspath(self.path)).split('_')
        name = f"{bits[0]}_{bits[-1]}"
        self.obj['name'] = name
        
        if self.runinfo:
            self.obj['RunInfo'] = self.runinfo.data
            if self.runinfo.recipe:
                self.obj['run_setup'] = self.runinfo.recipe
        
        if self.runparameters:
            self.obj.update(self.runparameters.data)
            if self.runparameters.recipe:
                self.obj['run_setup'] = self.runparameters.recipe
        
        if self.samplesheet:
            self.obj['samplesheet_csv'] = self.samplesheet.data
        
        if self.lanebarcodes:
            self.obj['illumina'] = {}
            self.obj['illumina']['Demultiplex_Stats'] = {}
            self.obj['illumina']['Demultiplex_Stats']['Barcode_lane_statistics'] = self.lanebarcodes.sample_data
            self.obj['illumina']['Demultiplex_Stats']['Flowcell_stats'] = self.lanebarcodes.flowcell_data
            if self.lanes:
                self.obj['illumina']['Demultiplex_Stats']['Lanes_stats'] = self.lanes.sample_data
        
        if self.interop_data:
            self.obj['InterOp'] = self.interop_data.data

        if self.undet:
            self.obj['Undetermined'] = self.undet.result
        
        if self.time_cycles:
            time_cycles = []
            for cycle in self.time_cycles.cycles:
                for k,v in cycle.items():
                    cycle[k] = str(v)
            self.obj['time cycles'] = self.time_cycles.cycles

        if self.json_stats:
            self.obj['Json_Stats'] = self.json_stats.data

class DemuxSummaryParser(object):
    def __init__(self, path):
        if os.path.exists(path):
            self.path = path
            self.result = {}
            self.TOTAL = {}
            self.parse()
        else:
            raise os.error(f"DemuxSummary folder {path} cannot be found")

    def parse(self):
        #will only save the 50 more frequent indexes
        pattern = re.compile('DemuxSummaryF1L([0-9]).txt')
        for file in glob.glob(os.path.join(self.path, 'DemuxSummaryF1L?.txt')):
            lane_nb = pattern.search(file).group(1)
            self.result[lane_nb] = OrderedDict()
            self.TOTAL[lane_nb] = 0
            with open(file, 'r') as f:
                undeterminePart = False
                for line in f:
                    if not undeterminePart:
                        if "### Columns:" in line:
                            undeterminePart = True
                    else:
                        #it means I am readng the index_Sequence  Hit_Count
                        components = line.rstrip().split('\t')
                        if len(self.result[lane_nb].keys())< 50:
                            self.result[lane_nb][components[0]] = int(components[1])
                        self.TOTAL[lane_nb] += int(components[1])


class LaneBarcodeParser(object):
    def __init__(self, path ):
        if os.path.exists(path):
            self.path = path
            self.parse()
        else:
            raise os.error(f" laneBarcode.html cannot be found at {path}")

    def parse(self):
        self.sample_data = []
        self.flowcell_data = {}
        with open(self.path, 'r') as htmlfile:
            bsoup = BeautifulSoup(htmlfile, 'html.parser')
            flowcell_table = bsoup.find_all('table')[1]
            lane_table = bsoup.find_all('table')[2]

            
            keys = []
            values = []
            for th in flowcell_table.find_all('th'):
                keys.append(th.text)
            for td in flowcell_table.find_all('td'):
                values.append(td.text)

            self.flowcell_data = dict(zip(keys, values))

            keys = []
            rows = lane_table.find_all('tr')
            for row in rows[0:]:
                if len(row.find_all('th')):
                    #this is the header row
                    for th in row.find_all('th'):
                        key = th.text.replace('<br/>', ' ').replace('&gt;', '>')
                        keys.append(key)
                elif len(row.find_all('td')):
                    values = []
                    for td in row.find_all('td'):
                        values.append(td.text)

                    d = dict(zip(keys,values))
                    self.sample_data.append(d)


class SampleSheetParser(object):
    """Parses  Samplesheets, with their fake csv format.

    Should be instancied with the samplesheet path as an argument.

    .header : a dict containing the info located under the [Header] section
    .settings : a dict containing the data from the [Settings] section
    .reads : a list of the values in the [Reads] section
    .data : a list of the values under the [Data] section. These values are stored in a dict format
    .datafields : a list of field names for the data section"""
    def __init__(self, path):
        self.log = logging.getLogger(__name__)
        if os.path.exists(path):
            self.parse(path)
        else:
            raise os.error(f"Sample sheet cannot be found at {path}")

    def parse(self, path):
        flag = None
        header = {}
        reads = []
        settings = {}
        csvlines = []
        data = []
        flag = 'data' #in case of HiSeq samplesheet only data section is present
        separator = ","
        with open(path, 'r') as csvfile:
            # Ignore empty lines (for instance the Illumina Experiment Manager
            # generates sample sheets with empty lines
            lines = filter(None, (line.rstrip() for line in csvfile))
            # Now parse the file
            for line in lines:
                if '[Header]' in line:
                    flag = 'HEADER'
                elif '[Reads]' in line:
                    flag = 'READS'
                elif '[Settings]' in line:
                    flag = 'SETTINGS'
                elif '[Data]' in line:
                    flag = 'data'
                else:
                    tokens = line.split(separator)
                    if flag == 'HEADER':
                        if len(tokens) < 2:
                            self.log.error("file {} does not seem has a correct format.")
                            raise RuntimeError("Could not parse the samplesheet, "
                                               "the file does not seem to have a correct format.")
                        header[tokens[0]] = tokens[1] 
                    elif flag == 'READS':
                        reads.append(tokens[0])
                    elif flag == 'SETTINGS':
                        settings[tokens[0]] = tokens[1]
                    elif flag == 'data':
                        csvlines.append(line)
       
            reader = csv.DictReader(csvlines)
            for row in reader:
                linedict = {}
                for field in reader.fieldnames:
                    linedict[field] = row[field]
                data.append(linedict)           
            self.datafields = reader.fieldnames
            self.dfield_sid = self._get_pattern_datafield(r'sample_?id')
            self.dfield_snm = self._get_pattern_datafield(r'sample_?name')
            self.dfield_proj = self._get_pattern_datafield(r'.*?project')
            self.data = data
            self.settings = settings
            self.header = header
            self.reads = reads

    def _get_pattern_datafield(self, pattern):
        for fld in self.datafields:
            if re.search(pattern,fld,re.IGNORECASE):
                return fld
        return ''


class SampleSheetV2Parser(object):
    """Parses the V2 Samplesheets, with their fake csv format.
    Should be instantiated with the samplesheet path as an argument.

    .header : a dict containing the info located under the [Header] section
    .reads : a list of the values in the [Reads] section
    .sequencer_settings : a dict containing the data from the [Sequencing_Settings] section
    .convert_settings : a dict containing the data from the [BCLConvert_Settings] section
    .cloud_settings : a dict containing the data from the [Cloud_Settings] section
    .convert_data : a list of the values under the [BCLConvert_Data] section. These values are stored in a dict format
    .convert_datafields : a list of field names for the [BCLConvert_Data] section
    .cloud_data : a list of the values under the [Cloud_Data] section. These values are stored in a dict format
    .cloud_datafields : a list of field names for the [Cloud_Data] section
    .data : a combination of cloud_data and convert_data, for legacy compatibility
    """
    def __init__(self, path:str):
        self.log = logging.getLogger(__name__)
        if os.path.exists(path):
            self.parse(path)
        else:
            raise os.error(f"Sample sheet cannot be found at {path}")
        
    def parse(self, path:str):
        flag = None
        header = {}
        reads = []
        convert_settings = {}
        sequencer_settings = {}
        cloud_settings = {}
        convert_data_fields = []
        convert_data = []
        cloud_data_fields = []
        cloud_data = []

        with open(path, 'r') as csvfile:
            # Ignore empty lines (for instance the Illumina Experiment Manager
            # generates sample sheets with empty lines
            lines = filter(None, (line.rstrip() for line in csvfile))
            # If header, set flag an go for next line
            for line in lines:
                if '[Header]' in line:
                    flag = 'HEADER'
                    continue
                elif '[Reads]' in line:
                    flag = 'READS'
                    continue
                elif '[BCLConvert_Settings]' in line:
                    flag = 'BCLCONVERT_SETTINGS'
                    continue
                elif '[Sequencing_Settings]' in line:
                    flag = 'SEQUENCER_SETTINGS'
                    continue
                elif '[Cloud_Settings]' in line:
                    flag = 'CLOUD_SETTINGS'
                    continue
                elif '[BCLConvert_Data]' in line:
                    flag = 'CONVERT_DATA'
                    continue
                elif '[Cloud_Data]' in line:
                    flag = 'CLOUD_DATA'
                    continue
                
                # Attempt to parse csv line 
                try:
                    tokens = next(csv.reader([line]))
                except csv.Error as e:
                    raise SyntaxError(f"Could not parse line '{line}' in {path}: {e}")
                
                if len(tokens) < 2:
                    self.log.error(f"File {path} does not seem has a correct format.")
                    raise SyntaxError("Could not parse the samplesheet, "
                                            "the file does not seem to have a correct format.")

                # Act on line depending based on current flag
                if flag == 'HEADER':
                    header[tokens[0]] = tokens[1]
                elif flag == 'READS':
                    reads.append(tokens[0])
                elif flag == 'BCLCONVERT_SETTINGS':
                    convert_settings[tokens[0]] = tokens[1]
                elif flag == 'SEQUENCER_SETTINGS':
                    sequencer_settings[tokens[0]] = tokens[1]
                elif flag == 'CLOUD_SETTINGS':
                    cloud_settings[tokens[0]] = tokens[1]
                elif flag == 'CONVERT_DATA':
                    if not convert_data_fields:
                        convert_data_fields = tokens
                    else:
                        convert_data.append(dict(zip(convert_data_fields, tokens)))
                elif flag == 'CLOUD_DATA':
                    if not cloud_data_fields:
                        cloud_data_fields = tokens
                    else:
                        cloud_data.append(dict(zip(cloud_data_fields, tokens)))
            
            # Set instance variables 
            self.header = header
            self.reads = reads
            self.sequencer_settings = sequencer_settings
            self.convert_settings = convert_settings
            self.cloud_settings = cloud_settings
            self.convert_data = convert_data
            self.convert_datafields = convert_data_fields
            self.cloud_data = cloud_data 
            self.cloud_datafields = cloud_data_fields
            self.dfield_sid = self._get_pattern_datafield(r'sample_?id')
            self.dfield_proj = self._get_pattern_datafield(r'project.*?')
            self.data = [a | b for a, b in zip(cloud_data, convert_data)]

    def _get_pattern_datafield(self, pattern:str):
        for fld in self.cloud_datafields:
            if re.search(pattern, fld,re.IGNORECASE):
                return fld
        return ''


class RunInfoParser(object):
    """Parses  RunInfo.xml.
    Should be instancied with the file path as an argument.

    .data : a list of hand-picked values :
     -Run ID
     -Run Number
     -Instrument
     -Flowcell name
     -Run Date
     -Reads metadata
     -Flowcell layout
    """
    def __init__(self, path ):
        self.data = {}
        self.recipe = None
        self.path = path
        if os.path.exists(path):
            self.parse()
        else:
            raise os.error(f" run info cannot be found at {path}")

    def parse(self):
        data = {}
        tree = ET.parse(self.path)
        root = tree.getroot()
        run = root.find('Run')
        data['Id'] = run.get('Id')
        data['Number'] = run.get('Number')
        data['Instrument'] = run.find('Instrument').text
        data['Flowcell'] = run.find('Flowcell').text
        data['Date'] = run.find('Date').text
        data['Reads'] = []
        for read in run.find('Reads').findall('Read'):
            data['Reads'].append(read.attrib)
        layout = run.find('FlowcellLayout')
        data['FlowcellLayout'] = layout.attrib
        self.data = data
        self.recipe = make_run_recipe(self.data.get('Reads', {}))

    def get_read_configuration(self):
        """return a list of dicts containig the Read Configuration
            """
        readConfig = []
        try:
            readConfig = self.data['Reads']
            return sorted(readConfig, key = lambda r: int(r.get("Number", 0)))
        except IOError:
            raise RuntimeError('Reads section not present in RunInfo. Check the FC folder.')
        

class RunParametersParser(object):
    """Parses a runParameters.xml file.
       This is a much more general xml parser, it will build a dict from the xml data.
       Attributes might be replaced if children nodes have the same tag as the attributes
       This does not happen in the current xml file, but if you're planning to reuse this, it may be of interest.
    """

    def __init__(self, path ):
        self.data = {}
        self.recipe = None
        self.path = path
        if os.path.exists(path):
            self.parse()
        else:
            raise os.error(f"RunParameters file cannot be found at {path}")
        
    def parse(self):
        data = {}
        tree = ET.parse(self.path)
        root = tree.getroot()
        self.data = xml_to_dict(root)
        self.recipe = make_run_recipe(self.data.get('Setup', {}).get('Reads', {}).get('Read', {}))


def make_run_recipe(reads):
    """Based on either runParameters of RunInfo, gathers the information as to how many
    readings are done and their length, e.g. 2x150"""
    nb_reads = 0
    nb_indexed_reads = 0
    numCycles = 0
    for read in reads:
        nb_reads += 1
        if read['IsIndexedRead'] == 'Y':
            nb_indexed_reads += 1
        else:
            if numCycles and numCycles != read['NumCycles']:
                logging.warning("NumCycles in not coherent")
            else:
                numCycles = read['NumCycles']

    if reads:
        return f"{nb_reads - nb_indexed_reads}x{numCycles}"
    return None

def xml_to_dict(root):
    current = None

    children = list(root)
    if children:
        current = {}
        duplicates = {}
        for child in children:
            if len(root.findall(child.tag))>1:
                if child.tag not in duplicates:
                    duplicates[child.tag] = []
                lower = xml_to_dict(child)
                duplicates[child.tag].extend(lower.values())
                current.update(duplicates)
            else:
                lower = xml_to_dict(child)
                current.update(lower)
    if root.attrib:
        if current:
            if [x in current for x in root.attrib]:
                current.update(root.attrib)
            else:
                current.update({'attribs':root.attribs})
        else:
            current = root.attrib
    if root.text and root.text.strip() != "":
        if current:
            if 'text' not in current:
                current['text'] = root.text
            else:
                #you're really pushing here, pal
                current['xml_text'] = root.text
        else:
            current = root.text
    return {root.tag:current}


class CycleTimesParser(object):
    def __init__(self, path):
        if os.path.exists(path):
            self.path = path
            self.cycles = []
            self.parse()
        else:
            raise os.error(f"file {path} cannot be found")

    def parse(self):
        """
        parse CycleTimes.txt and return ordered list of cycles
            CycleTimes.txt contains records: <date> <time> <barcode> <cycle> <info>
            one cycle contains a few records (defined by <cycle>)
            parser goes over records and saves the first record of each cycle as start time
            and the last record of each cycle as end time
        """
        data = []
        date_format = '%m/%d/%Y-%H:%M:%S.%f'
        with open(self.path, 'r') as file:
            cycle_times = file.readlines()
            # if file is empty, return
            if not cycle_times:
                return

            # first line is header, don't read it
            for cycle_line in cycle_times[1:]:
                # split line into strings
                cycle_list = cycle_line.split()
                cycle_time_obj = {}
                # parse datetime
                cycle_time_obj['datetime'] = datetime.strptime(f"{cycle_list[0]}-{cycle_list[1]}", date_format)
                # parse cycle number
                cycle_time_obj['cycle'] = int(cycle_list[3])
                # add object in the list
                data.append(cycle_time_obj)


        # take the first record as current cycle
        current_cycle = {
            'cycle_number': data[0]['cycle'],
            'start': data[0]['datetime'],
            'end': data[0]['datetime']
        }
        # compare each record with current cycle (except the first one)
        for record in data[1:]:
            # if we are at the same cycle
            if record['cycle'] == current_cycle['cycle_number']:
                # override end of cycle with current record
                current_cycle['end'] = record['datetime']
            # if a new cycle starts
            else:
                # save previous cycle
                self.cycles.append(current_cycle)
                # initialize new current_cycle
                current_cycle = {
                    'cycle_number': record['cycle'],
                    'start': record['datetime'],
                    'end': record['datetime']
                }
        # the last records is not saved inside the loop
        if current_cycle not in self.cycles:
            self.cycles.append(current_cycle)


class StatsParser(object):
    def __init__(self, path):
        if os.path.exists(path):
            self.path = path
            self.cycles = []
            self.data = None
            self.parse()
        else:
            raise FileNotFoundError(f"File {path} cannot be found")

    def parse(self):
        with open(self.path) as data:
            self.data = json.load(data)


class InterOpParser(object):
    """
    Parse tile data via the .bin files in the InterOp directory. 
    Includes much of the data otherwise only visible via Sequence Analysis Viewer.
    Relies on the interop package provided by Illumina.
    """
    def __init__(self, path):
        if os.path.exists(path):
            self.path = path
            self.data = None
            self.parse()
        else:
            raise FileNotFoundError(f"Directory {path} cannot be found")
    
    def parse(self):
        all_data = pd.DataFrame(imaging(self.path))
        all_data = all_data.fillna('')
        data_dict = all_data.to_dict()

        seen_tiles = set()
        condensed_data = []
        for index, tile in data_dict.get('Tile', {}).items():
            if tile not in seen_tiles:
                tile_data = {}
                tile_data['Tile'] = int(tile)
                tile_data['% Pass Filter'] = data_dict.get('% Pass Filter', {}).get(index, '')
                tile_data['% Occupied'] = data_dict.get('% Occupied', {}).get(index, '')
                seen_tiles.add(tile)
                condensed_data.append(tile_data)
        
        self.data = condensed_data
