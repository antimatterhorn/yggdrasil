from docutils import nodes
from docutils.parsers.rst import Directive

class YggdrasilLogoDirective(Directive):
    has_content = False

    def run(self):
        literal = nodes.literal_block(
            text=(
                "               _             _ _ \n"
                " _ _ ___ ___ _| |___ ___ ___|_| |\n"
                "| | | . | . | . |  _| .'|_ -| | |\n"
                "|_  |_  |_  |___|_| |__,|___|_|_|\n"
                "|___|___|___|       v0.8.5       \n"
            )
        )
        literal['language'] = 'text'  # tell Sphinx it's plain text
        return [literal]

def setup(app):
    app.add_directive("yggdrasil-logo", YggdrasilLogoDirective)
