import sys,os

#
# String Manipulation
#

def CSS(n,m=0):
            out = ''
            for number in range(n-m+1):
                out += str(n-number) + (',' if number != (n-m) else '')
            return out

def find_substring(string, substring):
    substring_length = len(substring)    
    def recurse(locations_found, start):
        location = string.find(substring, start)
        if location != -1:
            return recurse(locations_found + [location], location+substring_length)
        else:
            return locations_found
    return recurse([], 0)

#
# Template Handleling 
#

def seek_and_destroy(template,data,permissive=False):
    tmp_arg_mrk = '####'
    com_char = '//'
    special_chars = []
    try:
        template_file = open(template,'r')
    except:
        print '\33[31mError\33[0m: Template file:',template,'not found'
        sys.exit()
    filled_template = ''
    has_empty = False
    for template_line in template_file:
        if template_line.find(com_char)==0:
            continue
        template_args = find_substring(template_line,tmp_arg_mrk)
        line = template_line
        if len(template_args)%2 != 0:
            print '\33[31mError\33[0m: Template key marker error in:',template_line,'\n'
            print 'The correct usage is:'+tmp_arg_mrk+'Some Template Argument'+tmp_arg_mrk
            sys.exit()
        for index in range(len(template_args)/2):
            template_key = template_line[(template_args[2*index]+len(tmp_arg_mrk)):template_args[2*index+1]]
            template_pattern = tmp_arg_mrk+template_key+tmp_arg_mrk
            try:
                line = line.replace(template_pattern,data[template_key])
            except:
                print '\33[31mError\33[0m: Template key:',template_key,'was not found in dictionary'
                sys.exit()        
            has_empty = (data[template_key] == '')
        for char in special_chars:
            line = line.replace(char,str('\\'+char))
        filled_template += line
    if has_empty:
        wrn_msg = '\33[95mWarning\33[0m: '
        wrn_msg += 'data structure: '+str(data)+' has an empty template argument'
        if not permissive:
            wrn_msg += ', skkiped'
            filled_template = ""
        else:
            wrn_msg += ', written with missing arguments'
        print wrn_msg
    template_file.close()
    return filled_template

def MakeDir(Path):
        line = 'mkdir -p '+Path
        if os.path.exists(Path):
            print "\33[31mError\33[0m: Directory",Path,"already exists, please remove it first"
            sys.exit()
        os.system(line)
    
def WriteFile(Path,Content):
    F = open(Path,"w+")
    F.write(Content)
    F.close()

def CopyFile(Orig,Dest):
    cmd = 'cp '+Orig+' '+Dest
    os.system(cmd)