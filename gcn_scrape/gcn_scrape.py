from bs4 import BeautifulSoup
import requests

page = requests.get('https://gcn.gsfc.nasa.gov/selected.html')
soup = BeautifulSoup(page.content, 'lxml')

table = soup.find_all('dl')[0]
for row in table.find_all('dt'):
	link = row.find_all('a', href = True)
	link = link[0].get('href')
	full_link = f'https://gcn.gsfc.nasa.gov/{link}'
	gcn = requests.get(full_link)
	gcnsoup = BeautifulSoup(gcn.content, 'lxml')
	gcn_text = gcnsoup.find('body')
	

