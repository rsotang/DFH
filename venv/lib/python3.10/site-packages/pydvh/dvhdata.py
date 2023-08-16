from enum import Enum
import numpy as np
import copy


class DVHType(Enum):
    DDVH = 1
    CDVH = 2


class DVHData:
    def __init__(self, structure: str,
                 dvh_type: DVHType,
                 lower_dose: float,
                 upper_dose: float,
                 bin_width: float,
                 dose_array,
                 dose_unit,
                 volume_array,
                 volume_unit,
                 total_volume=None,
                 min_dose=None,
                 max_dose=None,
                 mean_dose=None,
                 modal_dose=None,
                 median_dose=None,
                 dose_coverage=None,
                 sampling_coverage=None):
        """

        :type bin_width: float
        :type upper_dose: float
        :type lower_dose: float
        :type dvh_type: DVHType -> DDVH, CDVH, None
        """
        self.sampling_coverage = sampling_coverage
        self.dose_coverage = dose_coverage
        self.median_dose = median_dose
        self.modal_dose = modal_dose
        self.mean_dose = mean_dose
        self.max_dose = max_dose
        self.min_dose = min_dose
        self.volume_unit = volume_unit
        self.dose_unit = dose_unit
        self.volume_array = volume_array
        self.dose_array = dose_array
        self._structure = structure
        self.bin_width = bin_width
        self.upper_dose = upper_dose + np.finfo(float).eps
        self.lower_dose = lower_dose
        self.dvh_type = dvh_type
        self.bin_number = int((upper_dose - lower_dose) / bin_width)
        self.overflow = 0.0
        self.underflow = 0.0
        if total_volume is None:
            self._total_volume = 0
        else:
            self._total_volume = total_volume

        if dose_array is None:
            self._initialize_dose_array(lower_dose, upper_dose, self.bin_width)
        if volume_array is None:
            self._initialize_volume_array(self.bin_number)

    def _initialize_dose_array(self, lower_dose_limit, upper_dose_limit, bin_width):
        bin_number = int((upper_dose_limit - lower_dose_limit) / bin_width)
        self.dose_array = np.zeros(bin_number)
        self.dose_array[0] = lower_dose_limit + (bin_width / 2.0)
        for i in range(1, bin_number):
            self.dose_array[i] = self.dose_array[i - 1] + bin_width

    def _initialize_volume_array(self, bin_number):
        self.volume_array = np.zeros(bin_number)

    def __repr__(self) -> str:
        return f"{self.dvh_type} for {self._structure}, Bins={self.bin_number}, " \
               f"Bin_width={self.bin_width}, structure_volume={self._total_volume}"

    def _add_entry(self, dose, volume):
        if dose > self.upper_dose:
            self.overflow += volume
        if dose < self.lower_dose:
            self.lower_dose += volume

        # find bin position for the dose
        pos = self._find_position_for_dose(dose)

        # update the volume_array
        self.volume_array[pos] += volume

        # update total volume of structure
        self._total_volume += volume

    def _find_position_for_dose(self, dose):
        return int((dose - self.lower_dose) / self.bin_width)

    def curve(self, with_volume_normalization=False):
        """A generator that returns tuples of (dose, volume). The dose is the mean value of 'bin'."""
        if with_volume_normalization:
            return zip(self.dose_array, self.volume_array * 100.0 / self._total_volume)
        else:
            return zip(self.dose_array, self.volume_array)

    def to_cumulative(self, eclipse_complied=True):
        if self.dvh_type == DVHType.CDVH:
            return self
        dvh_copy = copy.deepcopy(self)
        dvh_copy._convert_to_cumulative(eclipse_complied)
        return dvh_copy

    def to_differential(self, eclipse_complied=False):
        if self.dvh_type == DVHType.DDVH:
            return self
        dvh_copy = copy.deepcopy(self)
        dvh_copy._convert_to_differential(eclipse_complied)
        return dvh_copy

    def to_absolute_volume_values(self):
        dvh_copy = copy.deepcopy(self)
        if self.volume_unit == '%':
            dvh_copy.volume_array = dvh_copy.volume_array * dvh_copy.total_volume / 100
            # TODO: This is a quick fix. A volume class must be introduced storing values and units.
            dvh_copy.volume_unit = 'cm3'
        return dvh_copy

    def _convert_to_cumulative(self, eclipse_complied=False):
        self.volume_array[-1] += self.overflow
        for i in reversed(range(1, len(self.volume_array))):
            self.volume_array[i - 1] += self.volume_array[i]
        self.underflow += self.volume_array[0]
        self.dvh_type = DVHType.CDVH

        # shift dose / # vt[1] - ddvh.bin_width / 2 to comply with Eclipse exported files
        if eclipse_complied:
            for i in range(0, len(self.dose_array)):
                self.dose_array[i] = self.dose_array[i] - (self.bin_width / 2)

    def _convert_to_differential(self, eclipse_complied=False):
        self.underflow -= self.volume_array[0]
        for i in range(0, len(self.volume_array) - 1):
            self.volume_array[i] -= self.volume_array[i + 1]
        self.volume_array[-2] -= self.overflow
        self.dvh_type = DVHType.DDVH
        # shift dose / # vt[1] - ddvh.bin_width / 2 to comply with Eclipse exported files
        if eclipse_complied:
            for i in range(0, len(self.dose_array)):
                self.dose_array[i] = self.dose_array[i] + (self.bin_width / 2)

    def _get_position_in_volume_array(self, volume, volume_units):
        volume_in_dvh_units = self._get_volume_in_dvh_units(volume, volume_units)
        pos = len(self.volume_array) - np.searchsorted([*reversed(self.volume_array)], volume_in_dvh_units)
        return pos

    def _get_volume_in_dvh_units(self, volume, volume_units):
        if volume_units == self.volume_unit:
            return volume

        if volume_units == '%':
            if (self.volume_unit == 'cm3') or (self.volume_unit == 'mm3'):
                volume_in_absolute = volume * 0.01 * self.total_volume
                return volume_in_absolute

        if (volume_units == 'cm3') or (volume_units == 'mm3'):
            if self.volume_unit == '%':
                volume_in_relative = volume * 0.01 / self.total_volume
                return volume_in_relative

        if volume_units == 'cm3':
            if self.volume_unit == 'mm3':
                volume_in_mm3 = volume * 0.1 * 0.1 * 0.1
                return volume_in_mm3

        if volume_units == 'mm3':
            if self.volume_unit == 'cm3':
                volume_in_cm3 = volume / (0.1 * 0.1 * 0.1)
                return volume_in_cm3

    def get_dose_at_volume(self, volume, volume_units="%"):
        if self.dvh_type != DVHType.CDVH:
            raise ValueError("Operation is available on Cumulative DVH")

        pos = self._get_position_in_volume_array(volume, volume_units)

        if pos == 0:
            return self.volume_array[pos]

        v1 = self.volume_array[pos]
        v2 = self.volume_array[pos - 1]

        # y = a*x + b
        a = (v1 - v2) / self.bin_width
        b = v1 - a * self.dose_array[pos]
        volume_in_dvh_units = self._get_volume_in_dvh_units(volume, volume_units)
        return (volume_in_dvh_units - b) / a

    def get_volume_at_dose(self, dose):
        if self.dvh_type != DVHType.CDVH:
            raise ValueError("Operation is available on Cumulative DVH")
        pos = np.searchsorted(self.dose_array, dose)
        if pos == 0:
            return self.volume_array[pos]

        if pos == len(self.dose_array):
            return self.volume_array[pos-1]

        v1 = self.volume_array[pos]
        v2 = self.volume_array[pos - 1]

        # y = a*x + b
        a = (v1 - v2) / self.bin_width
        b = v1 - a * self.dose_array[pos]
        return a * dose + b

    @property
    def total_volume(self):
        if self._total_volume is None:
            if self.dvh_type == DVHType.DDVH:
                self._total_volume = np.sum(self.volume_array) + self.underflow + self.overflow

            if self.dvh_type == DVHType.CDVH:
                self._total_volume = self.volume_array[0] + self.underflow

        return self._total_volume

    @property
    def structure_name(self):
        return self._structure

    @classmethod
    def from_dose_matrix(cls, structure,
                         dvh_type,
                         voxel_doses,
                         dose_unit,
                         voxel_volume,
                         volume_unit,
                         upper_dose_limit,
                         lower_dose_limit=0,
                         bin_width=0.1):
        ddvh = cls(structure=structure,
                   dvh_type=DVHType.DDVH,
                   lower_dose=lower_dose_limit,
                   upper_dose=upper_dose_limit,
                   bin_width=bin_width,
                   dose_array=None,
                   volume_array=None,
                   dose_unit=dose_unit,
                   volume_unit=volume_unit
                   )

        for dose in voxel_doses:
            ddvh._add_entry(dose, voxel_volume)

        if dvh_type == DVHType.CDVH:
            return ddvh.to_cumulative(eclipse_complied=True)

        return ddvh

    @classmethod
    def from_DDVH_fixed_bin_width(cls,
                                  structure,
                                  dvh_type,
                                  dose_array: np.array,
                                  volume_array: np.array,
                                  dose_unit=None,
                                  volume_unit=None,
                                  total_volume=None,
                                  min_dose=None,
                                  max_dose=None,
                                  mean_dose=None,
                                  modal_dose=None,
                                  median_dose=None,
                                  dose_coverage=None,
                                  sampling_coverage=None):
        # Guard checks
        cls._validate_array_or_raise_exception(dose_array, "dose_array must be an array of float numbers.")
        cls._validate_array_or_raise_exception(volume_array, "volume_array must be an array of float numbers.")
        if len(dose_array) != len(volume_array):
            raise ValueError('dose_array and volume_array must have the same length.')

        dose_array = np.array(dose_array)
        volume_array = np.array(volume_array)
        bin_width = np.round(np.mean(np.diff(dose_array)), 5)

        if (total_volume is None) and (dvh_type == DVHType.DDVH):
            total_volume = np.sum(volume_array)

        return cls(structure=structure,
                   dvh_type=dvh_type,
                   lower_dose=0.0,
                   upper_dose=np.max(dose_array) + (bin_width / 2.0),
                   bin_width=bin_width,
                   dose_array=dose_array,
                   volume_array=volume_array,
                   dose_unit=dose_unit,
                   volume_unit=volume_unit,
                   total_volume=total_volume,
                   min_dose=min_dose,
                   max_dose=max_dose,
                   mean_dose=mean_dose,
                   modal_dose=modal_dose,
                   median_dose=median_dose,
                   dose_coverage=dose_coverage,
                   sampling_coverage=sampling_coverage)

    @classmethod
    def _validate_array_or_raise_exception(cls, varray, message):
        if not isinstance(varray, (np.ndarray, list)):
            raise ValueError(message)
