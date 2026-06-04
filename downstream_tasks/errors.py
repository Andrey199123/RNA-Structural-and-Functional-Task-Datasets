from __future__ import annotations


class DownstreamTasksError(Exception):
    """Base error for downstream-tasks tooling."""


class DownloadError(DownstreamTasksError):
    pass


class ExternalToolError(DownstreamTasksError):
    pass


class ValidationError(DownstreamTasksError):
    pass

