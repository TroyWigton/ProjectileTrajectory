#!/usr/bin/env python3
import sys
import re


def filter_math(content):
    # Regex for block math $$ ... $$ -> \f[ ... \f]
    # We use non-greedy match .*? and allow newlines with re.DOTALL
    # But wait, $$ is usually on its own line for blocks.
    # Doxygen \f[ must be used.

    # Replace block math first
    # Pattern: $$ (content) $$
    # Replacement: \f[\1\f]

    # We need to be careful not to match inside code blocks?
    # For a simple filter, we might assume nice formatting.

    # Handle Block Math $$ ... $$
    # Using a function for replacement to strip $$
    pattern_block = re.compile(r'\$\$(.*?)\$\$', re.DOTALL)
    content = pattern_block.sub(lambda m: r'\f[' + m.group(1) + r'\f]', content)

    # Handle Inline Math $ ... $
    # Negative lookbehind/ahead to avoid matching inside typical variable names if they used $var (but markdown usually doesn't)
    # But simplest is just $...$
    # We must exclude the case we just replaced (which now has \f[) — wait, no, the previous regex consumed the $$s

    # Pattern: $ (content) $
    # Replacement: \f$\1\f$
    # We should avoid matching across multiple lines for inline math usually.
    pattern_inline = re.compile(r'(?<!\\)\$(.*?)(?<!\\)\$')
    content = pattern_inline.sub(lambda m: r'\f$' + m.group(1) + r'\f$', content)

    return content


if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                content = f.read()
                # Apply filter only to markdown files
                if filename.lower().endswith('.md'):
                    sys.stdout.write(filter_math(content))
                else:
                    sys.stdout.write(content)
        except Exception:
            # If for some reason we fail, output nothing or empty
            pass
    else:
        # Read from stdin if no file (Doxygen sometimes pipes) // check doxy docs, usually it calls with filename
        pass
