import sys
import os
import inspect
import time
import requests
# sudo pip3 install requests

###################################################################
### download file class with sessions
###################################################################

class sessionDownloader():
	#=============================================================================================

	def __init__(self):
		self.session = requests.Session()

	def download_file(self,url,out_path,force=False,remake_age=3650,replace_empty=True):
		if os.path.exists(out_path):
			if force:
				print("Forced - Deleting ",out_path)
				os.remove(out_path)

			elif (time.time() - os.path.getctime(out_path))/60/60/24 > remake_age:
				print("Outdated - Deleting ",out_path,(time.time() - os.path.getctime(out_path))/60/60/24)
				os.remove(out_path)

			elif (replace_empty and os.path.getsize(out_path) == 0):
				os.remove(out_path)

		if not os.path.exists(out_path):
			try:
				resp = self.session.post(url)

				status = resp.status_code

				if status != 200:
					return {"status":"Error","status_code":resp.status_code,"details":resp}

				content = resp.text
				open(out_path,"w").write(content)
				return {"status":"Success"}
			except Exception as e:
				return {"status":"Error","error_type":str(e)}

if __name__ == "__main__":

	pdb_id = "2AST"
	pdb_file_path = os.path.dirname(inspect.stack()[0][1])

	url = "http://files.rcsb.org/download/" + pdb_id + ".pdb"
	out_path = os.path.abspath(os.path.join(pdb_file_path, pdb_id + ".pdb"))

	sessionDownloaderObj = sessionDownloader()
	response = sessionDownloaderObj.download_file(url,out_path)

	print(open(out_path).read())
	print(response)


    # parser.add_argument('--pdb_id', help="pdb_id to use", default='2r7e')
    # def download_pdb_file(pdb_id):
    # # Download pdb file
    # return pdbl.retrieve_pdb_file(
    #     pdb_code=pdb_id, pdir=LOCAL_STORAGE, file_format="pdb")