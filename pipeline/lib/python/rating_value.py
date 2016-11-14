import re

class RatingValue(object):
    ratval_re = re.compile(r"Rating name \(version\): (?P<name>.*?) " \
                            r"\(v(?P<version>\d+?)\)\n" \
                            r"Description: (?P<descr>.*?)\n" \
                            r"Value: (?P<value>.*?)$", \
                            flags=re.MULTILINE | re.DOTALL)

    def __init__(self, name, version, description, value):
        self.name = name
        self.version = version
        self.description = description
        self.value = value
        if self.value is None:
            self.valuestr = "%s" % self.value
        else:
            self.valuestr = "%.12g" % self.value

    def __str__(self):
        text  = "Rating name (version): %s (v%d)\n" % (self.name, self.version)
        text += "Description: %s\n" % self.description
        text += "Value: %s" % self.valuestr
        return text

    def get_short_string(self):
        return "%s (v%d): %s" % (self.name, self.version, self.valuestr)


def parse_string(string):
    """Parse a string for rating values and return a list
        of RatingValue objects.

        Input:
            string: The string to parse.

        Output:
            ratvals: A list of RatingValue objects gleaned from the string.
    """
    matches = RatingValue.ratval_re.finditer(string)
    ratvals = []
    for match in matches:
        grps = match.groupdict()
        if grps['value'] == "None":
            value = None
        else:
            value = float(grps['value'])
        ratvals.append(RatingValue(grps['name'], int(grps['version']), \
                                    grps['descr'], value))
    return ratvals


def read_file(ratfilenm):
    """Read a rating file, and return the rating values parsed from it.
    
        Input:
            ratfilenm: Name of the *.rat file.

        Output:
            ratvals: A list of RatingValue objects parsed from the file.
    """
    f = open(ratfilenm, 'r')
    contents = f.read()
    f.close()
    return parse_string(contents)
