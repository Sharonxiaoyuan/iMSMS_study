import sys

f1 = sys.argv[1]
f2 = sys.argv[2]

diff_thresh = 0.00002
max_diff = 0
with open(f1) as f1:
    with open(f2) as f2:
        while True:
            l1 = f1.readline()
            l2 = f2.readline()
            if l1 == "":
                break

            if l1 != l2:
                ss1 = l1[:-1].split("\t")
                ss2 = l2[:-1].split("\t")
                if len(ss1) != len(ss2):
                    print("Files differ in number of columns on these lines")
                    print(l1, l2)
                    break
                for i in range(len(ss1)):
                    not_warned = True
                    if ss1[i] != ss2[i]:
                        float_diff = abs(float(ss1[i]) - float(ss2[i]))
                        max_diff = max(max_diff, float_diff)
                        if float_diff > diff_thresh and not_warned:
                            print("Files differ significantly in value")
                            print(l1, l2)
                            print(i, ss1[i], ss2[i], "Float Diff:", float_diff)
                            not_warned = False

print("Maximum Difference:", max_diff)
