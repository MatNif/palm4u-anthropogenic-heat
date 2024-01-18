import os
import textwrap
import re
import yaml
import termcolor
import jinja2


project_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
docs_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

file_link_url = f'file://{project_dir}'

def add_palm_macros(env):
    env.macro(include_palm_namelist)
    env.macro(include_palm_logging_ids)
    env.macro(include_palm_output_quantities)
    env.macro(link_palm_repo_file)


def print_message_to_terminal(message_string, loglevel):
    message_wrapped = textwrap.fill(message_string, width=100, initial_indent=' ' * 9 + '-  ', subsequent_indent=' ' * 12)
    if loglevel == 'info':
        message = termcolor.colored('INFO:', 'green') + message_wrapped[5:]
    elif loglevel == 'warning':
        message = termcolor.colored('WARNING:', 'yellow') + message_wrapped[8:]
    elif loglevel == 'error':
        message = termcolor.colored('ERROR:', 'red') + message_wrapped[6:]
    else:
        raise ValueError(f'Unknown loglevel: "{loglevel}"')
    print(message)


def print_yaml_parameter_warning_to_terminal(path, parameter, text):
    print_message_to_terminal(
        message_string='reading yaml database "{}" with warning at parameter "{}". {}'.format(
            termcolor.colored(path, 'magenta'),
            termcolor.colored(parameter, 'cyan'),
            text,
        ),
        loglevel='warning'
    )


def print_yaml_parameter_error_to_terminal(path, parameter, text):
    print_message_to_terminal(
        message_string='reading yaml database "{}" fails at parameter "{}". {}'.format(
            termcolor.colored(path, 'magenta'),
            termcolor.colored(parameter, 'cyan'),
            text,
        ),
        loglevel='error'
    )


def validate_yaml_namelist_database(data_path, content_dict):
    """
    <parameter-name>
      category: ''
      type: ''
      shape: ()  # optional
      default:
        value:  # optional
        depends_on:  # optional
        value_of:  # optional
      si-unit:  # optional
      description:
        short: ''
        long: ''  # optional
      allowed_values:  # optional
        - value: 1.0
          description: ''
    """
    valid = True
    for parameter, parameter_dict in content_dict.items():

        if 'category' in parameter_dict and not 'categories' in parameter_dict:
            parameter_dict['categories'] = parameter_dict['category']
            #print_yaml_parameter_warning_to_terminal(
            #    data_path,
            #    parameter,
            #    'Deprecated field "category" is used. Please use field "categories" instead',
            #)
        if not 'categories' in parameter_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The field "categories" is mandatory but missing',
            )
        if isinstance(parameter_dict['categories'], str):
            parameter_dict['categories'] = [parameter_dict['categories']]
        if not isinstance(parameter_dict['categories'], list):
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The mandatory field "categories" must have a value of type "list"',
            )
            valid = False
        if not all(isinstance(v, str) for v in parameter_dict['categories']):
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The mandatory field "categories:" must have a value of type "list" '
                'which must contain items of type "str"',
            )
            valid = False

        if not 'type' in parameter_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The field "type" is mandatory but missing',
            )
            valid = False
        if not isinstance(parameter_dict['type'], str):
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The mandatory field "type" must have a value of type "str"',
            )
            valid = False
        if not bool(re.match('^C(?:\*\d+)?$|^I$|^L$|^R$|^D$', parameter_dict['type'])):
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The mandatory field "type" must have a value that starts with either C, C*<int>, I, L, R or D',
            )
            valid = False
        parameter_type_class = None
        if parameter_dict['type'].startswith('C'):
            parameter_type_class = str
        if parameter_dict['type'].startswith('I'):
            parameter_type_class = int
        if parameter_dict['type'].startswith('L'):
            parameter_type_class = bool
        if parameter_dict['type'].startswith('R'):
            parameter_type_class = float

        if 'shape' in parameter_dict:
            if not isinstance(parameter_dict['shape'], str):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The optional field "shape" must have a value of type "tuple". '
                    'An example for a shape of a 2D array would be "(10,50)"',
                )
                valid = False
            else:
                if not bool(re.match('^\((\d+,)*(\d+)\)$', parameter_dict['shape'])):
                    print_yaml_parameter_error_to_terminal(
                        data_path,
                        parameter,
                        'The optional field "shape" must have a value of type "tuple" which must contain items '
                        'of type "int". An example for a shape of a 2D array would be "(10,50)"',
                    )
                    valid = False
                else:
                    parameter_dict['shape'] = tuple([int(v) for v in parameter_dict['shape'][1:-1].split(',')])
                    if not all(isinstance(v, int) for v in parameter_dict['shape']):
                        print_yaml_parameter_error_to_terminal(
                            data_path,
                            parameter,
                            'The optional field "shape" must have a value of type "tuple" '
                            'which must contain items of type "int"',
                        )
                        valid = False

        if not 'default' in parameter_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The field "default" is mandatory but missing',
            )
            valid = False
        if not isinstance(parameter_dict['default'], dict):
            parameter_dict['default'] = dict(
                value=parameter_dict['default'],
            )
        if not any(map(lambda k: k in parameter_dict['default'].keys(), ['value', 'value_of', 'depends_on'])):
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The mandatory field "default" must contain a key '
                'that is either "value", "value_of" or "depends_on"',
            )
            valid = False

        if not isinstance(parameter_dict['default'], dict):
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The mandatory field "default" must have a value of type "dict"',
            )
            valid = False
        if 'value' in parameter_dict['default']:
            if parameter_dict['default']['value'] is not None:
                if parameter_type_class == bool:
                    if parameter_dict['default']['value'] in ['.T.', '.TRUE.']:
                        parameter_dict['default']['value'] = True
                    elif parameter_dict['default']['value'] in ['.F.', '.FALSE.']:
                        parameter_dict['default']['value'] = False
                    else:
                        print_yaml_parameter_error_to_terminal(
                            data_path,
                            parameter,
                            'The mandatory field "default.value" must have a value that represents '
                            'a Fortran LOGICAL and can be either .T., .TRUE., .F. or .FALSE.',
                        )
                        valid = False
                if 'shape' in parameter_dict and isinstance(parameter_dict['default']['value'], list):
                    if not all(isinstance(v, parameter_type_class) for v in parameter_dict['default']['value']):
                        print_yaml_parameter_error_to_terminal(
                            data_path,
                            parameter,
                            'The mandatory field "default.value" must have a value of type "list" '
                            'which must contain items of type "{}"'.format(str(parameter_type_class.__name__)),
                        )
                        valid = False
                else:
                    if not isinstance(parameter_dict['default']['value'], parameter_type_class):
                        print_yaml_parameter_error_to_terminal(
                            data_path,
                            parameter,
                            'The mandatory field "default.value" must have a value '
                            'of type "{}"'.format(str(parameter_type_class.__name__)),
                        )
                        valid = False
                if parameter_type_class == bool:
                    parameter_dict['default']['value'] = '.TRUE.' if parameter_dict['default']['value'] else '.FALSE.'
        elif 'value_of' in parameter_dict['default']:
            parameter_dict['default']['value'] = 'Value of {}'.format(parameter_dict['default']['value_of'])
        elif 'depends_on' in parameter_dict['default']:
            parameter_dict['default']['value'] = 'Depends on {}'.format(parameter_dict['default']['depends_on'])

        if 'si-unit' in parameter_dict:
            if not isinstance(parameter_dict['si-unit'], str):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The optional field "si-unit" must have a value of type "str"',
                )
                valid = False

        if not 'description' in parameter_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The field "description" is mandatory but missing',
            )
            valid = False
        if not isinstance(parameter_dict['description'], dict):
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The mandatory field "description" must have a value of type "dict"',
            )
            valid = False

        if not 'short' in parameter_dict['description']:
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The field "description.short" is mandatory but missing',
            )
            valid = False
        if not isinstance(parameter_dict['description']['short'], str):
            print_yaml_parameter_error_to_terminal(
                data_path,
                parameter,
                'The mandatory field "description.short" must have a value of type "str"',
            )
            valid = False

        if 'long' in parameter_dict['description']:
            if not isinstance(parameter_dict['description']['long'], str):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The optional field "description.long" must have a value of type "str"',
                )
                valid = False

        if 'allowed_values' in parameter_dict:
            if not isinstance(parameter_dict['allowed_values'], list):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The optional field "allowed_values" must have a value of type "list"',
                )
                valid = False
            if not all(isinstance(v, dict) for v in parameter_dict['allowed_values']):
                print_yaml_parameter_error_to_terminal(
                    data_path,
                    parameter,
                    'The optional field "allowed_values" must have a value of type "list" '
                    'which must contain items of type "dict"',
                )
                valid = False

            for allowed_value in parameter_dict['allowed_values']:
                if not 'value' in allowed_value:
                    print_yaml_parameter_error_to_terminal(
                        data_path,
                        parameter,
                        'The field "allowed_values[].value" is mandatory but missing',
                    )
                    valid = False
                if 'shape' in parameter_dict and isinstance(allowed_value['value'], list):
                    if not all(isinstance(v, parameter_type_class) for v in allowed_value['value']):
                        print_yaml_parameter_error_to_terminal(
                            data_path,
                            parameter,
                            'The mandatory field "allowed_values[].value" has a value of type "list" '
                            'which must contain items of type "{}"'.format(str(parameter_type_class.__name__)),
                        )
                        valid = False
                else:
                    if not isinstance(allowed_value['value'], parameter_type_class):
                        print_yaml_parameter_error_to_terminal(
                            data_path,
                            parameter,
                            'The mandatory field "allowed_values[].value" must have a value '
                            'of type "{}"'.format(str(parameter_type_class.__name__)),
                        )
                        valid = False

                if not 'description' in allowed_value:
                    print_yaml_parameter_error_to_terminal(
                        data_path,
                        parameter,
                        'The field "allowed_values[].description" is mandatory but missing',
                    )
                    valid = False
                if not isinstance(allowed_value['description'], str):
                    print_yaml_parameter_error_to_terminal(
                        data_path,
                        parameter,
                        'The mandatory field "allowed_values[].description" must have a value of type "str"',
                    )
                    valid = False
    return valid


def render_namelist_to_markdown_as_table(namelist, content_dict, link_table=True, link_path=''):
    format_str_head = '| {} | {} | {} |\n'
    format_str_body = '| {} | *{}* | {} |\n'
    output_str = format_str_head.format('Parameter', 'Default','Description')
    output_str += format_str_head.format('-', '-', '-', '-', '-')
    for parameter, parameter_dict in content_dict.items():
        output_str += format_str_body.format(
            '[{2}]({0}#{1}--{2})'.format(link_path, namelist, parameter) if link_table else parameter,
            parameter_dict['default']['value'] if parameter_dict['default']['value'] is not None else 'undefined',
            parameter_dict['description']['short'],
        )

    output_str += '\n\n'
    output_str += '*[I]: Integer\n'
    output_str += '*[I1]: Integer 8-Bit\n'
    output_str += '*[I2]: Integer 16-Bit\n'
    output_str += '*[I4]: Integer 32-Bit\n'
    output_str += '*[I8]: Integer 64-Bit\n'
    output_str += '*[R]: Real\n'
    output_str += '*[R1]: Real 8-Bit\n'
    output_str += '*[R2]: Real 16-Bit\n'
    output_str += '*[R4]: Real 32-Bit\n'
    output_str += '*[R8]: Real 64-Bit\n'
    output_str += '*[C]: Character\n'
    output_str += '*[L]: Logical\n'
    output_str += '*[D]: Derived Data-Type\n'
    output_str += '*[undefined]: This parameter has no usable default value and probably needs to be set by the user!\n'
    return output_str


def render_namelist_to_markdown(namelist, content_dict, heading_level=3):
    list_format_str = ': __{}:__ {}\n'
    output_str = '<br>\n'
    for parameter, parameter_dict in content_dict.items():
        output_str += '#' * heading_level + parameter + ' {#' + '{0}--{1}'.format(namelist, parameter) + '}\n\n'
        output_str += list_format_str.format(
            'Fortran Type', '{} {}'.format(
                parameter_dict['type'],
                '({})'.format(','.join([str(x) for x in parameter_dict['shape']])) if 'shape' in parameter_dict else '',
            )
        )
        output_str += list_format_str.format('Default', '*' + str(parameter_dict['default']['value']) + '*' if
        parameter_dict['default']['value'] is not None else '*undefined*')
        if 'si-unit' in parameter_dict:  # and re.search('^[I,R].*', parameter_dict['type']):
            output_str += list_format_str.format('SI-Unit', parameter_dict['si-unit'])
        output_str += '\n{}\n\n'.format(
            textwrap.indent(parameter_dict['description']['short'], ' ' * 4),
        )
        if 'long' in parameter_dict['description']:
            if parameter_dict['description']['long']:
                output_str += '{}\n\n'.format(
                    textwrap.indent(parameter_dict['description']['long'], ' ' * 4),
                )
        if 'allowed_values' in parameter_dict:
            output_str += '\n{}\n\n'.format(
                textwrap.indent('Currently {} choices are available:'.format(len(parameter_dict['allowed_values'])),
                                ' ' * 4),
            )
            for allowed_value in parameter_dict['allowed_values']:
                output_str += '    - *{}*\n\n{}\n\n'.format(
                    allowed_value['value'],
                    textwrap.indent(allowed_value['description'], ' ' * 8),
                )
        output_str += '\n<br>\n'
    output_str += '\n\n'
    output_str += '*[I]: Integer\n'
    output_str += '*[I1]: Integer 8-Bit\n'
    output_str += '*[I2]: Integer 16-Bit\n'
    output_str += '*[I4]: Integer 32-Bit\n'
    output_str += '*[I8]: Integer 64-Bit\n'
    output_str += '*[R]: Real\n'
    output_str += '*[R1]: Real 8-Bit\n'
    output_str += '*[R2]: Real 16-Bit\n'
    output_str += '*[R4]: Real 32-Bit\n'
    output_str += '*[R8]: Real 64-Bit\n'
    output_str += '*[C]: Character\n'
    output_str += '*[L]: Logical\n'
    output_str += '*[D]: Derived Data-Type\n'
    output_str += '*[undefined]: This parameter has no usable default value and probably needs to be set by the user!\n'
    return output_str


def include_palm_namelist(data_path, categories=['all'], as_table=False, link_table=True, link_path='', heading_level=3):
    call_string = 'include_palm_namelist(\'{}\', categories={}, as_table={}, link_table={}, link_path={}, heading_level={})'.format(
        data_path, categories, as_table, link_table, link_path, heading_level,
    )

    namelist = os.path.splitext(os.path.basename(data_path))[0]

    if isinstance(categories, str):
        categories = [categories]
    assert isinstance(categories, list)

    try:
        with open(os.path.join(docs_dir, 'content/data/namelists', data_path + '.yml')) as f:
            content_str = f.read()
    except OSError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to open YAML database file. See server terminal output for details.\n'

    try:
        j2env = jinja2.Environment()
        j2env.globals['link_palm_repo_file'] = link_palm_repo_file
        template = j2env.from_string(content_str)
        content_str = template.render()
    except Exception as e:
        print_message_to_terminal(
            message_string='processing yaml database with jinja2 "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to processing YAML database file with jinja2. See server terminal output for details.\n'

    try:
        content_dict = yaml.load(content_str, Loader=yaml.FullLoader)
    except yaml.parser.ParserError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to parse YAML database file. See server terminal output for details.\n'

    # valdate the yaml database content
    try:
        valid = validate_yaml_namelist_database(data_path, content_dict)
        if not valid:
            raise KeyError('namelist database validation failed')
    except KeyError:
        return '!!! warning\n    '+call_string+' failed! Error in YAML database layout. See server terminal for details.\n'

    # filter by category if required
    if 'all' not in categories:
        filtered_content_dict = dict()
        for parameter, parameter_dict in content_dict.items():
            if any(map(lambda c: c in parameter_dict['categories'], categories)):
                filtered_content_dict[parameter] = parameter_dict
        content_dict = filtered_content_dict

    if as_table:
        output_str = render_namelist_to_markdown_as_table(namelist, content_dict, link_table=link_table, link_path=link_path)
    else:
        output_str = render_namelist_to_markdown(namelist, content_dict, heading_level=heading_level)
    return output_str


def validate_yaml_logging_id_database(data_path, content_dict):
    """
    <logging_id>
      loglevel: ''
      message: ''
      description: ''
    """
    valid = True
    for logging_id, logging_id_dict in content_dict.items():

        if not 'loglevel' in logging_id_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                logging_id,
                'The field "loglevel" is mandatory but missing',
            )
            valid = False
        if not isinstance(logging_id_dict['loglevel'], str):
            print_yaml_parameter_error_to_terminal(
                data_path,
                logging_id,
                'The mandatory field "loglevel" must have a value of type "str"',
            )
            valid = False
        if not bool(re.match('^INFO$|^WARNING$|^ERROR$', logging_id_dict['loglevel'])):
            print_yaml_parameter_error_to_terminal(
                data_path,
                logging_id,
                'The mandatory field "loglevel" must have a value of either INFO, WARNING or ERROR',
            )
            valid = False

        if not 'message' in logging_id_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                logging_id,
                'The field "message" is mandatory but missing',
            )
            valid = False
        if not isinstance(logging_id_dict['message'], str):
            print_yaml_parameter_error_to_terminal(
                data_path,
                logging_id,
                'The mandatory field "message" must have a value of type "str"',
            )
            valid = False

        if not 'description' in logging_id_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                logging_id,
                'The field "description" is mandatory but missing',
            )
            valid = False
        if not isinstance(logging_id_dict['description'], str):
            print_yaml_parameter_error_to_terminal(
                data_path,
                logging_id,
                'The mandatory field "description" must have a value of type "str"',
            )
            valid = False
    return valid


def render_logging_ids_to_markdown_as_table(palm_module, content_dict, link_table=True, link_path=''):
    format_str_head = '| {} | {} | {} |\n'
    format_str_body = '| {} | *{}* | {} |\n'
    output_str = format_str_head.format('Logging ID', 'LogLevel', 'Message')
    output_str += format_str_head.format('-', '-', '-')
    for logging_id, logging_id_dict in content_dict.items():
        output_str += format_str_body.format(
            '[{2}]({0}#{1}--{2})'.format(link_path, palm_module, logging_id) if link_table else logging_id,
            logging_id_dict['loglevel'],
            str(logging_id_dict['message']).rstrip(),
        )

    output_str += '\n'
    return output_str


def render_logging_ids_to_markdown(palm_module, content_dict, heading_level=3):
    list_format_str = ': __{}:__ {}\n'
    output_str = ''
    for logging_id, logging_id_dict in content_dict.items():
        output_str += '#' * heading_level + logging_id + ' {#' + '{0}'.format(logging_id) + '}\n\n'
        output_str += list_format_str.format(
            'LogLevel', '*{}*\n\n'.format(logging_id_dict['loglevel'])
        )
        output_str += list_format_str.format(
            'Message', '{}\n\n'.format(
                textwrap.indent(logging_id_dict['message'], ' ' * 4),
            )
        )
        output_str += list_format_str.format(
            'Description', '{}\n\n'.format(
                textwrap.indent(logging_id_dict['description'], ' ' * 4),
            )
        )
    output_str += '\n'
    return output_str


def include_palm_logging_ids(data_path, loglevel=['all'], as_table=False, link_table=True, link_path='', heading_level=4):
    call_string = 'include_palm_logging_ids(\'{}\', loglevel={}, as_table={}, link_table={}, link_path={}, heading_level={})'.format(
        data_path, loglevel, as_table, link_table, link_path, heading_level,
    )

    palm_module = os.path.splitext(os.path.basename(data_path))[0]

    if isinstance(loglevel, str):
        loglevel = [loglevel]
    assert isinstance(loglevel, list)

    try:
        with open(os.path.join(docs_dir, 'content/data/logging', data_path + '.yml')) as f:
            content_dict = yaml.load(f.read(), Loader=yaml.FullLoader)
    except yaml.parser.ParserError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to parse YAML database file. See server terminal output for details.\n'
    except yaml.scanner.ScannerError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". This could be caused by using a ":" inside an unquoted string. {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to parse YAML database file. See server terminal output for details.\n'

    # valdate the yaml database content
    try:
        valid = validate_yaml_logging_id_database(data_path, content_dict)
        if not valid:
            raise KeyError('logging_ids database validation failed')
    except KeyError:
        return '!!! warning\n    '+call_string+' failed! Error in YAML database layout. See server terminal for details.\n'

    # filter by loglevel if required
    if 'all' not in loglevel:
        filtered_content_dict = dict()
        for logging_id, logging_id_dict in content_dict.items():
            if any(map(lambda c: c in logging_id_dict['loglevel'], loglevel)):
                filtered_content_dict[logging_id] = logging_id_dict
        content_dict = filtered_content_dict

    if as_table:
        output_str = render_logging_ids_to_markdown_as_table(palm_module, content_dict, link_table=link_table, link_path=link_path)
    else:
        output_str = render_logging_ids_to_markdown(palm_module, content_dict, heading_level=heading_level)
    return output_str


def validate_yaml_output_quantities_database(data_path, content_dict):
    """
    <output_quantity>
      scope:
        - palm_core
        - land_surface_model
      type:
        - vertical profile
        - 2d-array
        - 3d-array
        - masked array
      si-unit: %
      quantity: coverage of the land surface with bare soil
      remarks: ''
    """
    valid = True
    for output_quantity, output_quantity_dict in content_dict.items():

        if not 'scope' in output_quantity_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The field "scope" is mandatory but missing',
            )
            valid = False
        if not isinstance(output_quantity_dict['scope'], list):
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "scope" must have a value of type "list"',
            )
            valid = False
        if not all(isinstance(v, str) for v in output_quantity_dict['scope']):
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "scope" must have a value of type "list" '
                'which must contain items of type "str"',
            )
            valid = False

        if not 'type' in output_quantity_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The field "type" is mandatory but missing',
            )
            valid = False
        if not isinstance(output_quantity_dict['type'], list):
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "type" must have a value of type "list"',
            )
            valid = False
        if not all(isinstance(v, str) for v in output_quantity_dict['type']):
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "type" must have a value of type "list" '
                'which must contain items of type "str"',
            )
            valid = False
        if not all(bool(re.match('^vertical profile$|^2d-array$|^3d-array$|^masked array$', v)) for v in output_quantity_dict['type']):
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "type" must have a value of either "vertical profile", "2d-array", "3d-array" or "masked array"',
            )
            valid = False

        if isinstance(output_quantity_dict['si-unit'], int):
            output_quantity_dict['si-unit'] = str(output_quantity_dict['si-unit'])
        if not 'si-unit' in output_quantity_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The field "si-unit" is mandatory but missing',
            )
            valid = False
        if not isinstance(output_quantity_dict['si-unit'], str) and output_quantity_dict['si-unit'] != 1:
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "si-unit" must have a value of type "str"',
            )
            valid = False

        if not 'description' in output_quantity_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The field "description" is mandatory but missing',
            )
            valid = False
        if not isinstance(output_quantity_dict['description'], str):
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "description" must have a value of type "str"',
            )
            valid = False

        if not 'remarks' in output_quantity_dict:
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The field "remarks" is mandatory but missing',
            )
            valid = False
        if not isinstance(output_quantity_dict['remarks'], str):
            print_yaml_parameter_error_to_terminal(
                data_path,
                output_quantity,
                'The mandatory field "remarks" must have a value of type "str"',
            )
            valid = False
    return valid


def render_output_quantities_to_markdown_as_table(content_dict, show_remarks=True):
    format_str_head = '| {} | {} | {} |\n'
    format_str_body = '| {} | *{}* | {} |\n'
    output_str = format_str_head.format('Name { .oq-name }', 'SI-Unit { .oq-unit }', 'Description { .oq-description }')
    output_str += format_str_head.format('-', '-', '-')
    for output_quantity, output_quantity_dict in content_dict.items():
        output_quantity_excaped = output_quantity.replace('>', '\>').replace('_', '\_').replace('*', '\*')
        description_excaped = output_quantity_dict['description'].replace('>', '\>').replace('_', '\_').replace('*', '\*')
        remarks_excaped = output_quantity_dict['remarks'].replace('>', '\>').replace('_', '\_').replace('*', '\*')
        output_str += format_str_body.format(
            '{}'.format(output_quantity_excaped),
            output_quantity_dict['si-unit'],
            '{} {}'.format(
                str(description_excaped).rstrip(),
                '<br><br>*'+str(remarks_excaped).rstrip()+'*' if output_quantity_dict['remarks'] and show_remarks else '',
            ),
        )

    output_str += '\n'
    return output_str


def include_palm_output_quantities(data_path, scopes=['all'], types=['all'], show_remarks=True):
    call_string = 'include_palm_logging_ids(\'{}\', scopes={}, types={})'.format(
        data_path, scopes, types
    )

    palm_module = os.path.splitext(os.path.basename(data_path))[0]

    if isinstance(scopes, str):
        scopes = [scopes]
    assert isinstance(scopes, list)

    try:
        with open(os.path.join(docs_dir, 'content/data/output_quantities', data_path + '.yml')) as f:
            content_dict = yaml.load(f.read(), Loader=yaml.FullLoader)
    except yaml.parser.ParserError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to parse YAML database file. See server terminal output for details.\n'
    except yaml.scanner.ScannerError as e:
        print_message_to_terminal(
            message_string='reading yaml database "{}". This could be caused by using a ":" inside an unquoted string. {}'.format(
                termcolor.colored(data_path, 'magenta'),
                str(e)
            ),
            loglevel='error'
        )
        return f'!!! warning\n    {call_string} failed! Unable to parse YAML database file. See server terminal output for details.\n'

    # valdate the yaml database content
    try:
        valid = validate_yaml_output_quantities_database(data_path, content_dict)
        if not valid:
            raise KeyError('output_quantities database validation failed')
    except KeyError:
        return '!!! warning\n    '+call_string+' failed! Error in YAML database layout. See server terminal for details.\n'

    # filter by scope if required
    if 'all' not in scopes:
        filtered_content_dict = dict()
        for logging_id, logging_id_dict in content_dict.items():
            if any(map(lambda c: c in logging_id_dict['scope'], scopes)):
                filtered_content_dict[logging_id] = logging_id_dict
        content_dict = filtered_content_dict

    # filter by type if required
    if 'all' not in types:
        filtered_content_dict = dict()
        for logging_id, logging_id_dict in content_dict.items():
            if any(map(lambda c: c in logging_id_dict['type'], types)):
                filtered_content_dict[logging_id] = logging_id_dict
        content_dict = filtered_content_dict

    output_str = render_output_quantities_to_markdown_as_table(content_dict, show_remarks=show_remarks)
    return output_str


def link_palm_repo_file(link_text, file_path, new_tab=True):
    output_str = f'[{link_text}]({file_link_url}/{file_path})' + '{ target=_blank }' if new_tab else ''
    return output_str
