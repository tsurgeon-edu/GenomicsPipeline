import os
import subprocess
from datetime import datetime


class log_command(object):
    def __init__(self, command, from_function, th, function_class):
        self.command = command
        self.from_function = from_function
        self.th = str(th)
        self.f_class = function_class
        self.log = ""
        self.system_command_send()

    def system_command_send(self):
        start = datetime.now()
        # CSV-ish prefix: class,function,threads,start_time,
        prefix = f"{self.f_class},{self.from_function},{self.th},{start},"

        try:
            # Run the command in a shell so pipelines like "bwa ... | samtools ..." work
            result = subprocess.run(self.command, shell=True, text=True)

            end = datetime.now()
            if result.returncode == 0:
                log_line = prefix + f"{end},success," + self.command + "\n"
                self.write_logs(log_line)
                return True

            # Failure: log + stderr + raise
            log_line = prefix + f"{end},failed (exit_code={result.returncode})," + self.command + "\n"
            self.write_logs(log_line)

            if result.stderr:
                self.write_logs("STDERR:\n" + result.stderr + "\n")

            raise RuntimeError(
                f"{self.from_function} failed with exit code {result.returncode}\nCommand: {self.command}"
            )

        except Exception as e:
            end = datetime.now()
            log_line = prefix + f"{end},failed (exception),{self.command}\n"
            self.write_logs(log_line)
            # Re-raise so the pipeline stops and you see the real error
            raise

    def write_logs(self, log):
        with open("log_file.txt", "a") as file:
            file.write(log)
