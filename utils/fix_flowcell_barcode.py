
import GNomEx
g = GNomEx.GNomExConnection()
connection = g.GnConnect()
c = connection.cursor()
c.execute("update flowcell set barcode = 'C0C34ACXX' where barcode = 'COC34ACXX'")
x=connection.commit()
print "Commit returned %s" % x
