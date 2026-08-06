"""User-facing exceptions for :mod:`nucleic_builder`."""


class BuilderError(Exception):
    """Base class for expected input, conversion, and output errors."""


class InputValidationError(BuilderError):
    """The input is outside the selected canonical RNA or DNA contract."""


class ConversionError(BuilderError):
    """The pinned upstream converter could not produce a valid model."""


class OutputValidationError(BuilderError):
    """Generated molecule files are internally inconsistent."""


class OutputExistsError(BuilderError):
    """An output path already exists and overwrite was not requested."""
