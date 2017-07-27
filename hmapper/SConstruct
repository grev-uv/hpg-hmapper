
import os

# Initialize the environment with path variables, CFLAGS, and so on
bioinfo_path = '#lib/bioinfo-libs'
commons_path = '#lib/common-libs'

vars = Variables('buildvars.py')

compiler = ARGUMENTS.get('compiler', 'gcc')

env = Environment(tools = ['default', 'packaging'],
      		      CC = compiler,
                  variables = vars,
                  CFLAGS = '-std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp -msse4.2',
                  CPPPATH = ['#', '#src', '#include', bioinfo_path, commons_path, "%s/commons/argtable" % commons_path, "%s/commons/config" % commons_path ],
                  LIBPATH = ['#libs/common-libs/', commons_path],
                  LIBS = ['m', 'z', 'gsl', 'gslcblas'],
                  LINKFLAGS = ['-fopenmp'])
                  

if int(ARGUMENTS.get('debug', '0')) == 1:
    debug = 1
    env['CFLAGS'] += ' -DDEBUG -O0 -g'
else:
    debug = 0
    env['CFLAGS'] += ' -DDEBUG -O3 -g'

env['objects'] = []


# Targets

SConscript(['%s/SConscript' % bioinfo_path,
            '%s/SConscript' % commons_path
            ], exports = ['env', 'debug', 'compiler'])

env.Program('#bin/hpg-hmapper',
             source = [Glob('src/*.c'),
	               Glob('src/build-index/*.c'),
	               Glob('src/dna/*.c'),
	               Glob('src/rna/*.c'),
	               Glob('src/bs/*.c'),
                       "%s/libcommon.a" % commons_path,
                       "%s/libbioinfo.a" % bioinfo_path
                      ]
           )
