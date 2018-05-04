import requests
import xmltodict

def get_wufoo_textfile(url, zip_file=False):
    if zip_file == True:
        data = requests.get(url)
        return data
    else:
        data = requests.get(url)
        data.encoding = 'utf-8'
        data = data.text
        return data

def BBF_to_seq(URL):
    biobrick = xmltodict.parse(get_wufoo_textfile('http://parts.igem.org/cgi/xml/part.cgi?part=' + URL))
    part_type = biobrick["rsbpml"]["part_list"]["part"]["part_type"]
    print(part_type)
    seq = biobrick["rsbpml"]["part_list"]["part"]["sequences"]["seq_data"]
    return (seq,part_type)

