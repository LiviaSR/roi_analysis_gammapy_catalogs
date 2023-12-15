import sys
from io import StringIO
from datetime import datetime

class redirect_stdout(object):

    def __init__(self, write_func, **kwargs):
        """Init and set vars"""
        self.write_func = write_func
        self.kwargs = kwargs

    def __call__(self, function):
        """Call the function"""

        def wrapped_function(*args, **kwargs):

            try:
                sys.stdout = StringIO()
                function_to_exec = function(*args, **kwargs)
                out = sys.stdout.getvalue()
                self.write_func(out, **self.kwargs)

            finally:
                sys.stdout.close()
                sys.stdout = sys.__stdout__

            return function_to_exec

        return wrapped_function
    

def to_log_file(text, **kwargs):
    date = datetime.now().strftime("%y-%d-%m %H:%M:%S")
    file = kwargs["file"]
    file.write(f"{date}: \n{text}\n")
        