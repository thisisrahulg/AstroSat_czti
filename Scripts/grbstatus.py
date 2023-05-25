from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from bs4 import BeautifulSoup
from selenium.webdriver.chrome.options import Options
import sys

date = sys.argv[1]
time = sys.argv[2]
ra = sys.argv[3]
dec = sys.argv[4]

dt = f'{date} {time}'
rd = f'{ra}d {dec}d'

chrome_options = Options()
chrome_options.add_argument('--headless')
driverpath="/usr/bin/chromedriver"
service = Service(driverpath)

driver = webdriver.Chrome(service=service,options=chrome_options)
driver.get("http://astrosat-ssc.iucaa.in:8080/ASIMOV/ASIMOV.jsp")

datetime = driver.find_elements(By.ID,'DateTime')
datetime[0].send_keys(dt)

radec = driver.find_elements(By.ID,'coordinates')
radec[0].send_keys(rd)

astrosat = driver.find_elements(By.ID,'SatName_Check')
astrosat[0].click()

submit = driver.find_elements(By.ID,'SubmitButton')
submit[0].click()

page = driver.page_source

table = driver.find_elements(By.TAG_NAME,'table')
table_content = table[1].get_attribute('outerHTML')
soup = BeautifulSoup(table_content, 'html.parser')
column = soup.find_all('th')
ans = soup.find_all('td')
print('---------------------------------')
print(f'{column[1].text}:  {ans[1].text}')
print(f'{column[2].text}:  {ans[2].text}')
print('---------------------------------')
print(f'{column[3].text}:  {ans[3].text}')
print('---------------------------------')
driver.quit()
