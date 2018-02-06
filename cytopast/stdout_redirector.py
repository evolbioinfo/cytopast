import sys

import os


class HideOutput(object):
    """
    A context manager that block stdout for its scope, usage:

    with HideOutput():
        os.system('ls -l')
    """

    # One important thing to note is that a process gets three different file descriptors from the OS when launched:
    # stdin: 0
    # stdout: 1
    # stderr: 2

    def __init__(self, *args, **kw):
        sys.stdout.flush()
        self._origstdout = sys.stdout

        # Duplicate stdout to a different file descriptor number
        self._oldstdout_fno = os.dup(sys.stdout.fileno())

        # /dev/null is used just to discard what is being printed
        self._devnull = os.open(os.devnull, os.O_WRONLY)

    def __enter__(self):
        # Duplicate stdout (file descriptor 1) to a different file descriptor number
        self._newstdout = os.dup(1)

        # Duplicate the file descriptor for /dev/null and overwrite the value for stdout (file descriptor 1)
        os.dup2(self._devnull, 1)

        # Close devnull after duplication (no longer needed)
        os.close(self._devnull)

        # Use the original stdout to still be able to print to stdout within python
        sys.stdout = os.fdopen(self._newstdout, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._origstdout
        sys.stdout.flush()

        # Duplicate the file descriptor for old stdout and overwrite the value for stdout (file descriptor 1)
        os.dup2(self._oldstdout_fno, 1)
