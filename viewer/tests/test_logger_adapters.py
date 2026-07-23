"""Tests for ``viewer.logger_adapters.TaskLoggerAdapter`` (issue #991).

The adapter prefixes every log line with a compact task marker built from the
``extra`` dict. Each optional key (task, marker, target, tas, username) adds a
fragment; a blank username becomes ``u(|anon|)``. These are pure string
transforms - no logging output is captured, only ``process()`` is exercised.
"""
import logging

from viewer.logger_adapters import TaskLoggerAdapter


def _adapter(extra):
    return TaskLoggerAdapter(logging.getLogger("test"), extra)


def test_process_truncates_task_id():
    """Only the first 8 chars of the task UUID are rendered."""
    adapter = _adapter({"task": "54e3ea08-383a-4d96-9fc9-51948c8130b6"})
    msg, _ = adapter.process("hello", {})
    assert msg == "T.54e3ea08 hello"


def test_process_without_task_uses_placeholder():
    adapter = _adapter({"marker": "DOWNLOAD"})
    msg, _ = adapter.process("hello", {})
    assert msg == "T.000000 DOWNLOAD hello"


def test_process_all_fragments():
    adapter = _adapter(
        {
            "task": "abcdef0123456789",
            "marker": "LOAD",
            "target": "Mpro",
            "tas": "lb12345-1",
            "username": "alice",
        }
    )
    msg, _ = adapter.process("go", {})
    assert msg == "T.abcdef01 LOAD t(Mpro) tas(lb12345-1) u(alice) go"


def test_process_blank_username_is_anon():
    adapter = _adapter({"task": "abcdef0123456789", "username": ""})
    msg, _ = adapter.process("go", {})
    assert msg == "T.abcdef01 u(|anon|) go"


def test_process_passes_kwargs_through():
    adapter = _adapter({"task": "abcdef0123456789"})
    _, kwargs = adapter.process("go", {"stacklevel": 2})
    assert kwargs == {"stacklevel": 2}
