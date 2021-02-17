from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os
import pathlib
import shutil


class CMakeExtension(Extension):

    def __init__(self, name, sources=[]):
        super().__init__(name=name, sources=sources)


class BuildCMakeExt(build_ext):

    def run(self):

        for extension in self.extensions:
            if isinstance(extension, CMakeExtension):
                self.build_cmake(extension)

        self.extensions = list(filter(lambda e: not isinstance(e, CMakeExtension), self.extensions))
        super().run()

    def build_cmake(self, extension: Extension):

        self.announce("Preparing the build environment", level=3)

        workdir = os.getcwd()
        build_dir = pathlib.Path('build')
        os.makedirs(build_dir, exist_ok=True)

        extension_path = pathlib.Path(self.get_ext_fullpath(extension.name))
        install_path = extension_path.parent.absolute()
        extension_name = extension_path.parts[-1]
        os.makedirs(install_path, exist_ok=True)

        self.announce("configuring CMake project", level=3)
        os.chdir(build_dir)

        groops_env_var = os.getenv('GROOPS_SOURCE_DIR')
        if groops_env_var is None:
            groops_dir = pathlib.Path(os.path.join(os.path.expanduser('~'), 'groops', 'source'))
        else:
            groops_dir = pathlib.Path(groops_env_var)

        cmake_command = ['cmake', '..', '-DGROOPS_SOURCE_DIR={0}'.format(groops_dir), '-DEXTENSION_LIBRARY_NAME={0}'.format(extension_name)]
        if os.name == 'nt':
            cmake_command.extend(['-G', 'MinGW Makefiles'])
        self.spawn(cmake_command)

        self.announce("building binaries", level=3)
        self.spawn(["cmake", "--build", '.', "--target", extension_name, "--config", "Release"])
        shutil.copy(extension_name, install_path)
        os.chdir(workdir)


setup(
    name='groopsio',
    version='0.1',
    author='Andreas Kvas',
    description='A python package to read and write GROOPS files',
    long_description=open("README.md", 'r').read(),
    packages=['groopsio'],
    ext_modules=[CMakeExtension(name="groopsiobase")],
    license='GPL-3.0',
    cmdclass={'build_ext': BuildCMakeExt},
    install_requires=['cmake', 'numpy']
)
