import MySQLdb, config

class Database(object):
    def __init__(self, databasenm):
        self.database = config.DATABASES[databasenm]
        self.db = MySQLdb.connect(host=self.database["hostnm"],
                                  user=self.database["usernm"],
                                  db=self.database["dbnm"], 
                                  passwd=self.database["passwd"])
        self.cursor = self.db.cursor()

    def execute(self, query):
        try:
            self.cursor.execute(query)
        except Exception,e:
            print("ERROR: %i: %s"%(e.args[0], e.args[1]))
            self.db.rollback()
            self.db.close()
            exit(1)
    
    def commit(self):
        self.db.commit()

    def fetchone(self):
        return self.cursor.fetchone()

    def fetchall(self):
        return self.cursor.fetchall()
    
    def close(self):
        self.db.close()

def _strfmt(x):
    if type(x) == str and x != "NULL":
        return "\"%s\""%x
    else:
        return str(x)


