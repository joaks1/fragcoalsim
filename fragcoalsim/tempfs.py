#! /usr/bin/env python

import sys
import os
import logging
import tempfile

from fragcoalsim.stringutils import random_str
from fragcoalsim.errors import TempFSError

_LOG = logging.getLogger(__name__)

class TempFileSystem(object):
    """
    A temporary file system that protects against deleting any directories
    or files that are not created by an instance of this class.
    """
    def __init__(self, parent = None, prefix='fragcoalsim-temp-'):
        """
        Create TempFileSystem instance.

        `parent` must be an existing directory and will be where the base
        directory of the temp file system will be created. If not provided,
        tempfile.gettempdir() is used.
        """
        if not parent:
            parent = tempfile.gettempdir()
        if not os.path.exists(parent):
            raise TempFSError('{0!r} does not exist; cannot create '
                    'directory for temp files\n'.format(dest_dir))
        if not os.path.isdir(parent):
            raise TempFSError('{0!r} is not a directory; cannot create '
                    'directory for temp files\n'.format(dest_dir))
        self.dirs= set()
        self.files = set()
        self.parent = self._get_full_path(parent)
        self.prefix = prefix
        self.base_dir = self._make_dir(parent=self.parent, prefix=self.prefix)
        self.token_id = random_str()
        self.deleted = False

    def _get_full_path(self, path):
        return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))

    def _make_dir(self, parent, prefix):
        d = tempfile.mkdtemp(prefix=prefix, dir=parent)
        self._register_dir(d)
        return d

    def _register_dir(self, path):
        self.dirs.add(path)

    def _get_file(self, parent, prefix, register = True):
        prefix = '-'.join([self.token_id, prefix])
        file_descriptor, path = tempfile.mkstemp(prefix=prefix, dir=parent)
        os.close(file_descriptor)
        if register:
            self._register_file(path)
        return path
    
    def _register_file(self, path):
        self.files.add(path)

    def remove_file(self, path):
        if path in self.files:
            self.files.remove(path)
            os.remove(path)
        elif os.path.basename(path).startswith(self.token_id):
            os.remove(path)
        else:
            raise TempFSError('File {0!r} is not registered; '
                    'cannot remove'.format(path))

    def _check_parent(self, parent):
        full_parent = self._get_full_path(parent)
        if not os.path.exists(full_parent):
            raise TempFSError('parent {0!r} does not exist'.format(full_parent))
        if not os.path.isdir(full_parent):
            raise TempFSError('parent {0!r} is not a directory'.format(
                    full_parent))
        if not full_parent in self.dirs:
            raise TempFSError('unregistered parent: {0}'.format(full_parent))
        return full_parent

    def get_file_path(self, parent=None, prefix='temp', register = True):
        """
        Get temp file path within the temp directory.

        `parent` must exist within temp file base directory and must
        be registered with this TempFileSystem instance.

        if `register` is True, the temp file is registered (i.e., it is added
        to the instance's set of temp files).
        """
        if parent is None:
            parent = self.base_dir
        full_parent = self._check_parent(parent)
        return self._get_file(parent=full_parent, prefix=prefix, register = register)

    def create_subdir(self, parent=None, prefix='temp'):
        """
        Create a new subdirectory within the temp file system.

        `parent` must exist within temp file base directory and must
        be registered with this TempFileSystem instance.
        """
        if parent is None:
            parent = self.base_dir
        full_parent = self._check_parent(parent)
        return self._make_dir(parent=full_parent, prefix=prefix)
        
    def clear_dir(self, path):
        full_path = self._get_full_path(path)
        if not full_path in self.dirs:
            raise TempFSError('Temp directory {0!r} is not registered; '
                    'cannot clear'.format(full_path))
        for p in os.listdir(full_path):
            path = os.path.join(full_path, p)
            if os.path.isfile(path):
                self.remove_file(path)
            else:
                self.remove_dir(path)

    def remove_dir(self, path):
        full_path = self._get_full_path(path)
        self.clear_dir(full_path)
        self.dirs.remove(full_path)
        try:
            os.rmdir(full_path)
        except OSError as e:
            _LOG.warning('Could not remove temp directory {0}. Here are the '
                    'contents:\n{1}'.format(full_path,
                            '\n'.join(os.listdir(full_path))))
            pass
        if full_path == self.base_dir:
            self.deleted = True
                
    def purge(self):
        self.remove_dir(self.base_dir)

