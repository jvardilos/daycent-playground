with open("aggregations.txt", "r") as file:
    content = file.read()

listing = content.split("\n")

# get total from all som
aggregate = 0
for list in listing:
    if list != "":
        sepBySpace = list.split(" ")
        aggregate += float(sepBySpace[0])


print('total = {}'.format(aggregate))

valueList = []
percentList = []
for list in listing:
    if list != "":
        sepBySpace = list.split(" ")
        key = sepBySpace[len(sepBySpace)-1]
        value = float(sepBySpace[0])
        if value != 0.0:
            valueList.append(value)
            percentVal = value / aggregate * 100
            percentList.append("{} = {}".format(key, percentVal))

for value in valueList:
    print(value)


aggregatePercent = 0
for list in percentList:
    print(list)
    l = list.split(" = ")
    aggregatePercent += float(l[len(l)-1])

print('total percent = {}'.format(aggregatePercent))
