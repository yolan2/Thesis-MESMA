
import re

filename = r'c:\Users\yolan\OneDrive\Documenten\UGENT\Master\R_MESMA\fit_veg_mixture_mesma.R'
try:
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()
except UnicodeDecodeError:
    with open(filename, 'r', encoding='latin-1') as f:
        content = f.read()

# Remove comments (naive)
# content_no_comments = re.sub(r'#.*', '', content)
# Better comment removal (handling strings)
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
                out.append(' ') # Replace comment with space to keep line numbers
        elif in_string:
            if char == string_char:
                if i > 0 and text[i-1] == '\\':
                    pass # Escaped quote
                else:
                    in_string = False
            out.append(' ') # Replace string content with space
        else:
            if char == '#' and (i == 0 or text[i-1] != '\\'): # Naive check for escaped #
                in_comment = True
                out.append(' ')
            elif char == '"' or char == "'" or char == '`':  # Include backtick for R function names
                in_string = True
                string_char = char
                out.append(' ')
            else:
                out.append(char)
        i += 1
    return "".join(out)

clean_content = remove_comments_and_strings(content)

stack = []
lines = clean_content.split('\n')
matching = {'}': '{', ')': '(', ']': '['}

for i, line in enumerate(lines):
    for j, char in enumerate(line):
        if char in '{([':
            stack.append((i + 1, j + 1, char))
        elif char in '})]':
            if not stack:
                print(f'Extra closing {char} at line {i + 1}, col {j + 1}')
            else:
                line_num, col_num, opening = stack[-1]
                if opening == matching[char]:
                    stack.pop()
                else:
                    print(f'Mismatched delimiter: {opening} at line {line_num}, col {col_num} closed by {char} at line {i + 1}, col {j + 1}')

if stack:
    print(f'Unclosed delimiters: {len(stack)}')
    for item in stack[-5:]:
        print(f'Unclosed {item[2]} at line {item[0]}, col {item[1]}')
else:
    print('All delimiters balanced.')
