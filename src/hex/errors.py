import sys


def terminate_program(function_name: str, info_message: str) -> None:
    print(info_message)
    message = 'Critical error at "' + function_name + '"-function; program terminated'
    sys.exit(message)