from typing import Union, Tuple


class NameConversionError(Exception):
    def __init__(self, message, lipid):
        super(NameConversionError, self).__init__(message)

        self.message = message
        self.lipid = lipid


class LipidDataError(Exception):
    def __init__(self, message):
        super(LipidDataError, self).__init__(message)


class GroupDataError(Exception):
    def __init__(self, message):
        super(GroupDataError, self).__init__(message)


class FaSettingError(Exception):
    def __init__(self, message, error_type):
        super(FaSettingError, self).__init__(message)
        self.error_type = error_type
        self.message = message


class LipidSettingError(Exception):
    def __init__(self, message, error_type):
        super(LipidSettingError, self).__init__(message)
        self.error_type = error_type
        self.message = message


class MissingClassError(Exception):
    def __init__(self, message, class_, lipid, error_location):
        super(MissingClassError, self).__init__(message)
        self.message = message
        self.class_ = class_
        self.lipid = lipid
        self.location = error_location


class NotComputedError(Exception):
    def __init__(self, attribute_name: str,
                 function_name: Union[str, Tuple[str]] = None,
                 for_subset: bool = False):
        if function_name is None:
            super(NotComputedError, self).__init__(
                attribute_name
            )
        else:
            if for_subset:
                message = f"{attribute_name} has not been computed for {function_name}"
            else:
                message = f"{attribute_name} has not been computed yet,"\
                          f" please call {function_name} first!"
            super(NotComputedError, self).__init__(message)


class CorrelationError(Exception):
    def __init__(self, message, lipids):
        super(CorrelationError, self).__init__(message)

        self.lipids = lipids
        self.message = message


class PartialCorrelationError(Exception):
    def __init__(self, message, solver):
        super(PartialCorrelationError, self).__init__(message)

        self.solver = solver
        self.message = message


class SignificanceTestError(Exception):
    def __init__(self, message):
        super(SignificanceTestError, self).__init__(message)


class SignificanceTestNAError(Exception):
    def __init__(self, message):
        super(SignificanceTestNAError, self).__init__(message)
        self.message = message
