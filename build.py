import re
import argparse
from jinja2 import Environment, FileSystemLoader


THIS_DIR = '.'

user_struct = [
    ('outdir', 'char[256]', '.'),
    ('tmax', 'double', 1.0),
    ('cpi', 'double', 1.0),
    ('N', 'int[2]', [128, 1]),
    ('ng', 'int[2]', [0, 0]),
    ('domain', 'double[4]', [0, 1, 0, 1]),
    ('coordinates', 'char', 'p'),
    ('move_cells', 'int', 0),
    ('initial_data', 'char[256]', 'uniform'),
    ('initial_mesh', 'char[256]', 'planar'),
    ('boundary_con', 'char[256]', ''),
]


status_struct = [
    ('iteration', 'int', 0),
    ('checkpoint_number', 'int', 0),
    ('time_simulation', 'double', 0.0),
    ('time_step', 'double', 0.0),
    ('time_last_checkpoint', 'double', 0.0),
    ('kzps', 'double', 0.0),
]


j2_env = Environment(loader=FileSystemLoader(THIS_DIR),
                     trim_blocks=True)


def make_member(name, dtype, default):
    match = re.match('(.*)\[([0-9]+)\]', dtype)
    typ = match.groups(0)[0] if match else dtype
    arr = int(match.groups(0)[1]) if match else None
    mem = dict(name=name, dtype=typ)
    if not arr:
        mem['printf_arg'] = [name]
        mem['printf_fmt'] = {'double':'%4.3lf',
                             'int':'%d',
                             'char':'%c'}[typ]
    elif typ == 'char':
        mem['printf_arg'] = [name]
        mem['printf_fmt'] = '%s'
    elif typ == 'int':
        mem['printf_arg'] = [name + '[%d]'%n for n in range(arr)]
        mem['printf_fmt'] = '%d ' * arr
    elif typ == 'double':
        mem['printf_arg'] = [name + '[%d]'%n for n in range(arr)]
        mem['printf_fmt'] = '%4.3lf ' * arr
    mem['array'] = arr
    mem['default_value'] = default
    return mem



def render(infiles):
    template_env = dict(app_name='gusto',
                        app_use_hdf5=True,
                        user_struct=[make_member(*m) for m in user_struct],
                        status_struct=[make_member(*m) for m in status_struct])
    for infile in infiles:
        template = j2_env.get_template(infile)
        print template.render(**template_env)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infiles', nargs='+')
    args = parser.parse_args()
    render(args.infiles)
