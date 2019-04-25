import lxml.html

player_page = requests.get('http://www.naic.edu/~pfreire/GCpsr.html')
player_tree = lxml.html.fromstring(player_page.content)

print player_tree