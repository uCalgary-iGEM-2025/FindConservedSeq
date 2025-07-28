# ERROR CHECK IN CASE VARIABLES ARE EMPTY

import parse
import mutations

instances = parse.parse_nucleotides("test.fasta")
top_row = parse.construct_top_row(instances)
print(top_row)
for i in instances:
    print(i)
print()
new_instances, new_top_row = mutations.shift_mutations(instances)
print(new_top_row)
for i in new_instances:
    print(i)
