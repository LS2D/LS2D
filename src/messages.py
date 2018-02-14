_opts = {
   'blue'   : '\033[94m',
   'green'  : '\033[92m',
   'purple' : '\033[95m',
   'red'    : '\033[91m',
   'bf'     : '\033[1m',
   'ul'     : '\033[4m',
   'end'    : '\033[0m'
}

def header(message):
    """
    Format of print statements indicating new main routine
    """
    print('{}{}{}{}'.format(_opts['blue'], _opts['bf'], message, _opts['end']))

def message(message):
    """
    Format of print statements
    """
    print(' - {}'.format(message))

def warning(message):
    """
    Format of print warnings
    """
    print('{}{}WARNING:{} {}'.format(_opts['purple'], _opts['bf'], _opts['end'], message))

def error(message):
    """
    Format of print errors
    """
    print('{}{}ERROR:{} {}'.format(_opts['red'], _opts['bf'], _opts['end'], message))
