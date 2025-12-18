import re
fn = r'c:\Users\yolan\OneDrive\Documenten\UGENT\Master\R_MESMA\fit_veg_mixture_mesma.R'

def remove_comments_and_strings(text):
    out = []
    i = 0
    n = len(text)
    in_string = False
    string_char = ''
    in_comment = False
    while i < n:
        char = text[i]
        if in_comment:
            if char == '\n':
                in_comment = False
                out.append(char)
            else:
                out.append(' ')
        elif in_string:
            if char == string_char:
                if i > 0 and text[i-1] == '\\':
                    pass
                else:
                    in_string = False
            out.append(' ')
        else:
            if char == '#' and (i == 0 or text[i-1] != '\\'):
                in_comment = True
                out.append(' ')
            elif char in '"\'`':
                in_string = True
                string_char = char
                out.append(' ')
            else:
                out.append(char)
        i += 1
    return ''.join(out)

with open(fn,'r',encoding='utf-8') as f:
    s = f.read()
clean = remove_comments_and_strings(s)
lines = clean.split('\n')
stack=[]
for i,line in enumerate(lines, start=1):
    for j,ch in enumerate(line, start=1):
        if ch == '{':
            stack.append((i,j,ch))
            print(f'PUSH {{ at {i}:{j} (stack depth {len(stack)})')
        elif ch == '}':
            if stack and stack[-1][2] == '{':
                o = stack.pop()
                print(f'POP  }} at {i}:{j} matched with {o[0]}:{o[1]} (stack depth {len(stack)})')
            else:
                print(f'UNMATCHED }} at {i}:{j} (no matching {{ on stack)')

print('remaining stack size:', len(stack))
if stack:
    for item in stack[-20:]:
        print('UNclosed { at', item)
