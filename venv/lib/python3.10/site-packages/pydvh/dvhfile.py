from pathlib import Path
from io import open
import re
import numpy as np
from pydvh.dvhdata import DVHType, DVHData


class DVHFile:

    def __init__(self, filename, patient_name=None, patient_id=None, type: DVHType = None,
                 course_name=None,
                 plan_sum_name=None,
                 plan_name=None,
                 structure_names=None,
                 dvhs=None) -> None:
        super().__init__()
        self._dvhs = dvhs
        self._structure_names = structure_names
        self._plan_name = plan_name
        self._plan_sum_name = plan_sum_name
        self._course_name = course_name
        self._type = type
        self._patient_id = patient_id
        self._patient_name = patient_name
        self.filename = filename

    def patient_name(self):
        return self._patient_name

    def get_dvh_by_name(self, structure):
        dvhs = [*(x for x in self._dvhs if x._structure == structure)]
        if not dvhs:
            return None
        else:
            return dvhs[0]

    def structure_names(self):
        return self._structure_names

    def course_name(self):
        return self._course_name

    def type(self):
        return self._type

    def patient_id(self):
        return self._patient_id

    def plan_sum_name(self):
        return self._plan_sum_name

    def plan_name(self):
        return self._plan_name

    @classmethod
    def from_file_eclipse(cls, filename):
        """
        Parses VARIAN Eclipse DVH exported file. Eclipse v15.6 tested.
        :param filename: Eclipse DVH exported file (v15.6)
        :return: DVHFile object which contains the DVHData for all available structure.
        """
        dvhfilename = Path(filename)
        if not dvhfilename.is_file():
            return DVHFileNull(filename)
        dvhfile = cls.__open_dvhfile(filename)
        lines = [*dvhfile.readlines()]
        patient_name = cls._parse_patient_name(lines)
        patient_id = cls._parse_patient_id(lines)
        dvhfile_type = cls._parse_type(lines)
        course_name = cls._parse_course(lines)
        plan_sum_name = cls._parse_plan_sum(lines)
        if plan_sum_name is None:
            plan_name = cls._parse_plan(lines)
        else:
            plan_name = None
        structure_names = cls._parse_structure_names(lines)
        dvhs = cls._get_dvhs_eclipse(structure_names, lines, dvhfile_type)

        dvhfile.close()
        return DVHFile(filename=filename,
                       patient_name=patient_name,
                       patient_id=patient_id,
                       type=dvhfile_type,
                       course_name=course_name,
                       plan_sum_name=plan_sum_name,
                       plan_name=plan_name,
                       structure_names=structure_names,
                       dvhs=dvhs)

    @classmethod
    def __open_dvhfile(cls, filename):
        return open(filename, "r")

    @classmethod
    def _parse_patient_name(cls, lines):
        for line in lines:
            res = re.findall('Patient Name(.*): \\s*([^\n\r]*)', line)
            if len(res) > 0:
                return res[0][1]
        return "NULL"

    @classmethod
    def _parse_patient_id(cls, lines):
        for line in lines:
            res = re.findall('Patient ID(.*): \\s*([^\n\r]*)', line)
            if len(res) > 0:
                return res[0][1]
        return "NULL"

    @classmethod
    def _parse_type(cls, lines):
        for line in lines:
            res = re.findall('^Type(.*): \\s*([^\n\r]*)', line)
            if len(res) > 0:
                if "Differential" in res[0][1]:
                    return DVHType.DDVH
                if "Cumulative" in res[0][1]:
                    return DVHType.CDVH
        return None

    @classmethod
    def _parse_course(cls, lines):
        for line in lines:
            res = re.findall('Course(.*): \\s*([^\n\r]*)', line)
            if len(res) > 0:
                return res[0][1]
        return None

    @classmethod
    def _parse_plan_sum(cls, lines):
        for line in lines:
            res = re.findall('Plan sum(.*): \\s*([^\n\r]*)', line)
            if len(res) > 0:
                return res[0][1]
        return None

    @classmethod
    def _parse_plan(cls, lines):
        for line in lines:
            res = re.findall('Plan: \\s*([^\n\r]*)', line)
            if len(res) > 0:
                return res[0]
        return None

    @classmethod
    def _parse_structure_names(cls, lines):
        names = []
        for line in lines:
            res = re.findall('Structure: \\s*([^\n\r]*)', line)
            if len(res) > 0:
                names.append(res[0])
        return names

    @classmethod
    def _parse_get_line_number_structure(cls, structure, lines):
        names = []
        counter = 0
        for line in lines:
            res = re.findall('Structure: \\s*([^\n\r]*)', line)
            if len(res) > 0:
                if res[0] == structure:
                    return counter
            counter += 1
        return names

    @classmethod
    def _get_dvhs_eclipse(cls, structure_names, lines, dvh_type):
        partitions = cls._get_partitions(structure_names, lines)
        dvhs = []
        for partition in partitions:
            s = partition['structure']
            dvh = cls._get_dvh_eclipse(structure=s,
                                       partition=partition,
                                       lines=lines,
                                       dvh_type=dvh_type)
            dvhs.append(dvh)
        return dvhs

    @classmethod
    def _get_partitions(cls, structure_names, lines):
        partitions = []
        for structure in structure_names:
            start = cls._parse_get_line_number_structure(structure, lines)
            stop = start
            for i in np.arange(start + 1, len(lines)):
                stop = i - 1
                if "Structure:" in lines[i]:  # the next structure reached
                    break
            partitions.append({'structure': structure,
                               'start': start,
                               'stop': stop})
        return partitions

    @classmethod
    def _get_dvh_eclipse(cls, structure: str, partition, lines, dvh_type: DVHType) -> object:
        dose_array = []
        volume_array = []
        start = partition['start']
        stop = partition['stop']
        can_start_gather_dose_volume = False
        for i in range(start, stop):
            data = lines[i]
            if can_start_gather_dose_volume == True:
                cols = re.split('[\s]+', data)
                if len(cols) == 4:
                    dose_array.append(float(cols[1]))
                    volume_array.append(float(cols[2]))
                if len(cols) == 5:
                    dose_array.append(float(cols[1]))
                    volume_array.append(float(cols[3]))

            if dvh_type == DVHType.DDVH:
                if ("Dose [Gy] dVolume" in data) \
                        or ("Dose [Gy]   Relative dose [%]" in data):
                    can_start_gather_dose_volume = True
            if dvh_type == DVHType.CDVH:
                if ("Dose [Gy] Ratio of Total Structure Volume [%]" in data) or (
                        "Dose [Gy]   Relative dose [%]" in data):
                    can_start_gather_dose_volume = True

        dose_unit = cls._get_dose_unit(lines, partition)
        volume_unit = cls._get_volume_unit(lines, partition)
        dose_array = np.array(dose_array)
        volume_array = np.array(volume_array)

        bin_width = np.mean(np.diff(dose_array))
        structure_volume = None

        if dvh_type == DVHType.DDVH:
            volume_array = np.array(volume_array) * bin_width
        if dvh_type == DVHType.CDVH:
            volume_array = np.array(volume_array)
            structure_volume = cls._get_structure_volume(lines, partition)
        min_dose = cls._get_minimum_dose(lines, partition)
        max_dose = cls._get_maximum_dose(lines, partition)
        mean_dose = cls._get_mean_dose(lines, partition)
        modal_dose = cls._get_modal_dose(lines, partition)
        median_dose = cls._get_median_dose(lines, partition)
        dose_coverage = cls._get_dose_coverage(lines, partition)
        sampling_coverage = cls._get_sampling_coverage(lines, partition)
        return DVHData.from_DDVH_fixed_bin_width(structure=structure,
                                                 dvh_type=dvh_type,
                                                 dose_array=dose_array,
                                                 volume_array=volume_array,
                                                 dose_unit=dose_unit,
                                                 volume_unit=volume_unit,
                                                 total_volume=structure_volume,
                                                 min_dose=min_dose,
                                                 max_dose=max_dose,
                                                 mean_dose=mean_dose,
                                                 modal_dose=modal_dose,
                                                 median_dose=median_dose,
                                                 dose_coverage=dose_coverage,
                                                 sampling_coverage=sampling_coverage)

    @classmethod
    def _get_dose_unit(cls, lines, partition):
        dict = {'%': 0,
                'Gy': 0,
                'cGy': 0}
        start = partition['start']
        stop = partition['stop']
        for i in range(start, stop):
            if "Dose" in lines[i]:
                if "%" in lines[i]:
                    dict['%'] += 1
                if "Gy" in lines[i]:
                    dict['Gy'] += 1
                if "cGy" in lines[i]:
                    dict['cGy'] += 1

        max_k = max(dict, key=dict.get)
        return max_k

    @classmethod
    def _get_volume_unit(cls, lines, partition):
        dict = {'%': 0,
                'cm3': 0}
        start = partition['start']
        stop = partition['stop']
        for i in range(start, stop):
            if "Volume" in lines[i]:
                if "%" in lines[i]:
                    dict['%'] += 1
                if "cm" in lines[i]:
                    dict['cm3'] += 1

        max_k = max(dict, key=dict.get)
        return max_k

    @classmethod
    def _get_structure_volume(cls, lines, partition):
        start = partition['start']
        stop = partition['stop']
        for i in range(start, stop):
            res = re.findall('Volume(.*): \\s*([^\n\r]*)', lines[i])
            if len(res) > 0:
                return float(res[0][1])
        return None

    @classmethod
    def _get_minimum_dose(cls, lines, partition):
        start = partition['start']
        stop = partition['stop']
        for i in range(start, stop):
            res = re.findall('Min Dose(.*): \\s*([^\n\r]*)', lines[i])
            if len(res) > 0:
                return float(res[0][1])
        return None


    @classmethod
    def _get_maximum_dose(cls, lines, partition):
        start = partition['start']
        stop = partition['stop']
        for i in range(start, stop):
            res = re.findall('Max Dose(.*): \\s*([^\n\r]*)', lines[i])
            if len(res) > 0:
                return float(res[0][1])
        return None

    @classmethod
    def _get_mean_dose(cls, lines, partition):
        start = partition['start']
        stop = partition['stop']
        for i in range(start, stop):
            res = re.findall('Mean Dose(.*): \\s*([^\n\r]*)', lines[i])
            if len(res) > 0:
                return float(res[0][1])
        return None


    @classmethod
    def _get_modal_dose(cls, lines, partition):
        start = partition['start']
        stop = partition['stop']
        for i in range(start, stop):
            res = re.findall('Modal Dose(.*): \\s*([^\n\r]*)', lines[i])
            if len(res) > 0:
                return float(res[0][1])
        return None

    @classmethod
    def _get_median_dose(cls, lines, partition):
        start = partition['start']
        stop = partition['stop']
        for i in range(start, stop):
            res = re.findall('Median Dose(.*): \\s*([^\n\r]*)', lines[i])
            if len(res) > 0:
                return float(res[0][1])
        return None

    @classmethod
    def _get_dose_coverage(cls, lines, partition):
        start = partition['start']
        stop = partition['stop']
        for i in range(start, stop):
            res = re.findall('Dose Cover(.*): \\s*([^\n\r]*)', lines[i])
            if len(res) > 0:
                return float(res[0][1])
        return None

    @classmethod
    def _get_sampling_coverage(cls, lines, partition):
        start = partition['start']
        stop = partition['stop']
        for i in range(start, stop):
            res = re.findall('Sampling Cover(.*): \\s*([^\n\r]*)', lines[i])
            if len(res) > 0:
                return float(res[0][1])
        return None


class DVHFileNull(DVHFile):
    def __init__(self, filename) -> None:
        super().__init__(filename=filename, patient_name="NULL")
