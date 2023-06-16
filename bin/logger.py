import datetime

log_file = './log/error.log'

def log(message):
    with open(log_file, 'a') as f:
        current_time = datetime.datetime.now().strftime('%Y/%m/%d %H:%M')
        f.write(f'{current_time}\t{message}\n')