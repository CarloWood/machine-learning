directory /usr/src/Arch/tensorflow/tensorflow/src/tensorflow-2.14.0-opt-cuda-dbg

set pagination off

# skip all STL source files
define skipstl
python
# get all sources loadable by gdb
def GetSources():
    sources = []
    for line in gdb.execute('info sources',to_string=True).splitlines():
        if line.startswith("/"):
            sources += [source.strip() for source in line.split(",")]
    return sources

# skip files of which the (absolute) path begins with 'dir'
def SkipDir(dir):
    sources = GetSources()
    for source in sources:
        if source.startswith(dir):
            gdb.execute('skip file %s' % source, to_string=True)

# apply only for c++
if 'c++' in gdb.execute('show language', to_string=True):
    SkipDir("/usr")
    SkipDir("/opt")
end
end

define lvinit
    skipstl
    skip operator new(unsigned long)
end

define hookpost-run
    lvinit
end
