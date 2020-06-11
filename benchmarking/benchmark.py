from subprocess import check_output

template_file = "robot_template.txt"
output_file = "basic/base.vxa"
template_str = open(template_file, "r").read()

log_file = "data.csv"

LOOPS_PER_SIZE = 1

times = []

if __name__ == "__main__":
    for robot_len_pow in range(16):
        robot_len = 1<<robot_len_pow
        layer_str = "1"*3 *robot_len
        curr_robot = template_str % (str(robot_len), layer_str, layer_str, layer_str)
        with open(output_file, "w") as f:
            f.write(curr_robot)

        # after file is written, benchmark it
        for loop_id in range(LOOPS_PER_SIZE):
            print(str(robot_len*9), ",", sep="", end="", flush=True)
            time_line = [l for l in check_output("./voxcraft-sim -i basic".split()).decode("utf-8").split("\n") if "seconds" in l][0]
            took_idx = time_line.find("took") + 5
            seconds_idx = time_line.find("seconds.") - 1
            curr_time = float(time_line[took_idx:seconds_idx])
            times.append((robot_len*9, curr_time))
            print("%4.3f" %(curr_time))

    with open(log_file, "w") as f:
        for size, t in times:
            f.write("%d,%.3f\n"%(size, t))
    print(times)