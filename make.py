import os
import re
import sys

include_re = re.compile(r'\s*#\s*include\s*"(.*)"\s*')
main_re = re.compile(r'\s*int\s*main\s*\(.*\)\s*{?\s*')
makefile = open('Makefile', 'w')
include_dict = {}
main = 'main'
main_suffix = ''
del_command = 'rm'
if sys.platform == 'win32':
    main_suffix = '.exe'
    del_command = 'del'
test_file_re = re.compile(r'test\d*.cpp')

for root, dirs, filenames in os.walk("."):
    for filename in filenames:
        if test_file_re.match(filename):
            continue
        origin_name, suffix = os.path.splitext(os.path.join(root, filename))
        if suffix == ".cpp" or suffix == ".cc":
            include_dict[origin_name] = set()
            check_list = []
            check_list.append(os.path.join(root, filename))
            while check_list:
                fp = open(check_list[0])
                dir_name = os.path.dirname(check_list[0])
                for line in fp:
                    if main_re.match(line):
                        main = os.path.splitext(os.path.basename(check_list[0]))[0]
                    include_match = include_re.match(line)
                    if include_match:
                        header_name = include_match.group(1).strip()
                        if include_match.group(1) not in include_dict[origin_name]:
                            header_file = os.path.join(dir_name, header_name)
                            if os.path.isfile(header_file):
                                include_dict[origin_name].add(os.path.relpath(header_file, '.'))
                                check_list.append(header_file)
                            else:
                                print 'dir_name:%s' % dir_name
                                print 'header_name:%s' % header_name
                                print "%s isn't a file" % header_file
                fp.close()
                check_list.remove(check_list[0])
value_dict = {}
value_dict['cc'] = 'g++'
value_dict['exe'] = main + main_suffix
value_dict['std'] = 'c++11'
value_dict['objects'] = ' '.join(map(lambda x: x + '.o', include_dict.keys()))

choices = ['-O2', '-fopenmp']

for k, v in value_dict.items():
    makefile.write('%s = %s\n' % (k, v))
makefile.write('''
$(exe) : $(objects)
\t\t\t$(cc) %s -std=$(std) -o $(exe) $(objects)

''' % ' '.join(choices))

for k, v in include_dict.items():
    makefile.write('%s.o : %s.cpp %s\n' % (k, k, ' '.join(v)))
    makefile.write('\t$(cc) -std=$(std) %s -c %s.cpp -o %s\n' % (' '.join(choices), k, k + '.o'))

makefile.write('''
clean :
\t\t%s $(objects) $(exe)
''' % del_command)
