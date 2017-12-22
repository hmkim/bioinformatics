import string

term_list='''
18S ribosomal RNA
18S small subunit ribosomal RNA
25S ribosomal RNA
26S ribosomal RNA
28S large subunit ribosomal RNA
28S ribosomal RNA
28S ribosomal RNA, large subunit
5.8S ribosomal RNA
contains large subunit ribosomal RNA, internal transcribed spacer, and small subunit ribosomal RNA
contains small subunit ribosomal RNA, internal transcribed spacer, and large subunit ribosomal RNA
internal transcribed spacer
internal transcribed spacer 1
internal transcribed spacer 2
internal transcribed spacer; ITS
large subunit ribosomal RNA
small subunit ribosomal RNA
small subunit ribosomal RNA and internal transcribed spacer 1
'''

def main():
    new = term_list.replace('and', '\n')
    new = term_list.replace(';', '\n')
    print(new)

if __name__ == '__main__':
    main()    

