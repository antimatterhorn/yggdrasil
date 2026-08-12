from yggdrasil import *
from Calculators import TimeDilation
from Colors import BRIGHT_CYAN, BRIGHT_GREEN, BRIGHT_YELLOW, BRIGHT_RED, colorize

if __name__ == "__main__":
    commandLine = CommandLineArguments(travelTime=2, #years
                                        distance=4   #light-years
                                        )


    td = TimeDilation(travelTime,distance)

    print("To travel",colorize("%3.2f light-years"%distance,BRIGHT_YELLOW),"in",
        colorize("%3.2f years"%travelTime,BRIGHT_YELLOW),"(subective time),",
        "\nyou would need to travel ",colorize("%3.3f c"%td.velocity,BRIGHT_RED),
        "\nand",colorize("%3.3f years"%td.earthTime,BRIGHT_CYAN),
        "would pass on Earth")