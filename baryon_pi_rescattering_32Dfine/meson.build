project('qlat-cpp', 'cpp',
  version: '0.1',
  license: 'GPL-3.0-or-later',
  default_options: [
    'warning_level=0',
    'cpp_std=c++17',
    'libdir=lib',
    'optimization=2',
    'debug=false',
    ])

add_project_arguments('-fno-strict-aliasing', language: ['c', 'cpp'])
add_project_arguments('-O4', '-g3','-fno-strict-aliasing', language: ['c', 'cpp'])

cpp = meson.get_compiler('cpp')

omp = dependency('openmp')

# mpic = dependency('mpi', language: 'cpp')
# mpic = dependency('mpich', required : true)

cxx = run_command('bash', '-c', 'echo "$CXX"', check: true).stdout().strip()
mpicxx = run_command('bash', '-c', 'echo "$MPICXX"', check: true).stdout().strip()
mpic = dependency('mpi', language: 'cpp').as_system()


qlat = dependency('qlat')

deps = [ omp, mpic, qlat, ]

# deps = [ omp, qlat, ]

c = run_command('bash', '-c', 'ls "$MESON_SOURCE_ROOT"/*.cpp', check: true)
sources = c.stdout().strip().split('\n')

executable('qlat.x',
  sources,
  dependencies: deps,
  install: true,
  )

run_target('run',
  command: 'run.sh')
