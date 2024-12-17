class NomalyBaseException(Exception):
    """Base exception class for Nomaly application"""
    pass

class DatabaseConnectionError(NomalyBaseException):
    """Raised when database connection fails"""
    pass

class DataNotFoundError(NomalyBaseException):
    """Raised when requested data is not found"""
    pass

class GWASError(NomalyBaseException):
    """Raised when GWAS analysis fails"""
    pass

class PheWASError(NomalyBaseException):
    """Raised when PheWAS analysis fails"""
    pass

class ValidationError(NomalyBaseException):
    """Raised when input validation fails"""
    pass 