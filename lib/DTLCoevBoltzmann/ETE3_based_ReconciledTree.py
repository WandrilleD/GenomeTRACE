#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################################
##  Author:         Wandrille Duchemin  #
##  Created:        13-Jan-2017         #
##  Last modified:  16-Jan-2017         #
#########################################


import ete3 

EVENTTAGCORRESPONDANCE = {	"D" : "duplication",
						"S" : "speciation",
						"C" : "leaf",
						"SL": "speciationLoss",

						"Bo": "bifurcationOut",
						"So": "speciationOut",
						"Tb": "transferBack",
						"SoL": "speciationOutLoss"
						}

class RecEvent:
	def __init__(self, eventCode , species , ts = None , additionnalInfo = None):
		"""
		Takes:
			- eventCode (str) : a code indicating the recEvent event 
			- species (~) : a identifier for the specie the event takes place in
			- ts (int or None) [default= None] : the time slice the events happens at, if applicable
			- additionnalInfo (dict or None) [default= None] : keys are expected to be some property tag and values the associated information
		"""
		self.eventCode = eventCode
		self.species = species
		self.timeSlice = ts
		self.additionnalInfo = additionnalInfo

	def __str__(self):
		eventName = self.eventCode
		if EVENTTAGCORRESPONDANCE.has_key(self.eventCode):
			eventName = EVENTTAGCORRESPONDANCE[self.eventCode]

		L = [ str(self.eventCode), "spe=" + str(self.species) ]
		if not self.timeSlice is None:
			L.append("ts=" + str(self.timeSlice))
		for k,v in self.additionnalInfo.items():
			L.append(str(k) + "=" + str(v) )
		return " ".join(L)

	def nwkstr(self):
		""" tmp simplistic version """
		s = str(self.species)
		s += "."
		coev = self.additionnalInfo.get("coevent",False)
		if coev:
			if coev > 1:
				## co-TL case -> differentiate coTL , TcoL and coTcoL
				if coev == 2:
					s += self.eventCode[:-1] + "co" + self.eventCode[-1]
					return s 
				elif coev == 3:
					s += "co" + self.eventCode[:-1] + "co" + self.eventCode[-1]
					return s
			s += "co"
		s += str(self.eventCode)
		return s

	def isCoev(self):
		return self.additionnalInfo.get("coevent",False)

	def makeRecXMLstr(self , speciesNames):
		if self.eventCode == "N":
			return ""

		eventName = self.eventCode
		if EVENTTAGCORRESPONDANCE.has_key(self.eventCode):
			eventName = EVENTTAGCORRESPONDANCE[self.eventCode]

		spe = str(self.species)
		if spe in speciesNames.keys():
			spe = speciesNames[spe]

		S = "<" + str(eventName)

		if self.eventCode != "Bo":
			S += " "
			if self.eventCode == "Tb":
				S += "destinationSpecies="
			else:
				S += "speciesLocation=" 
			S += '"' + str(spe) + '"'

			if not self.timeSlice is None:
				S += " ts=" + '"' + str(self.timeSlice) + '"'

		S += ">"
		S += "</" + eventName + ">"
		
		return S


class ReconciledTree(ete3.TreeNode):
	def __init__(self):
		ete3.TreeNode.__init__(self)

		self.name = ""
		self.eventRecs = []

	def setName(self, name):
		""" NB : name can be any object with a __str__() fc (like an int for instance) """
		self.name = name

	def getEvents(self):
		return self.eventRecs

	def addEvent(self , e,append = True):
		"""
			Takes:
				- e (RecEvent) : the reconciliation event to add
				- append (bool) [default=True] : if True adds the event at the end of the list else, adds it at the beginning
		"""
		if append:
			self.eventRecs.append(e)
		else:
			self.eventRecs.insert(0,e)
		return

	def getTreeStrAux(self):

		L = []
		L.append("name : " + str(self.name) )
		L.append("events :")
		for e in self.eventRecs:
			L.append("  " + str(e))
		ChL = []
		for c in self.get_children():
			ChL += ["  " + s for s in c.getTreeStrAux()]
		if len(ChL)>0:
			L.append("children :")
			L += ChL
		return L



	def getTreeStr(self):
		return "\n".join(self.getTreeStrAux())

	def getTreeNewickAux(self, sep="|", topoOnly = False):
		s = ""

		ChL = []
		for c in self.get_children():
			ChL.append(c.getTreeNewickAux(sep,topoOnly))
		if len(ChL)>0:
			s += "("
			s += ",".join(ChL)
			s += ")"

		s += str(self.name)
		if not topoOnly:
			s += sep
			s += sep.join( [e.nwkstr() for e in self.eventRecs] )

		return s

	def getTreeNewick(self, sep="|" , topoOnly = False):
		return self.getTreeNewickAux(sep, topoOnly) + ";"

	def getTreeRecPhyloXMLAux(self , speciesNames ={}, topoOnly = False):

		L = []
		L.append("<clade>")
		L.append("  <name>" + str(self.name) + "</name>" )
		if not topoOnly:
			L.append("  <eventsRec>")
			for e in self.eventRecs:
				s = e.makeRecXMLstr(speciesNames)
				if s == "":
					continue
				L.append("    " +  s)
			L.append("  </eventsRec>")
		ChL = []
		for c in self.get_children():
			ChL += ["  " + s for s in c.getTreeRecPhyloXMLAux(speciesNames, topoOnly)]
		if len(ChL)>0:
			L += ChL
		L.append("</clade>")
		return L

	def getTreeRecPhyloXML(self , speciesNames ={}, topoOnly = False):
		Lines = ["<recGeneTree>"]
		Lines.append("  <phylogeny rooted=\"true\">")
		tmp = self.getTreeRecPhyloXMLAux( speciesNames , topoOnly )
		for l in tmp:
			Lines.append( "    " + l )
		Lines.append("  </phylogeny>")
		Lines.append( "</recGeneTree>" )

		return "\n".join(Lines)


	def countEvents(self, accountForCoEvent = False):

		devent = {}
		for e in self.eventRecs:
			code  = e.eventCode
			if accountForCoEvent:
				if e.isCoev():
					code  = "co" + code
			if not devent.has_key(code):
				devent[code] = 0
			devent[code] += 1 

		for c in self.get_children():
			tmp = c.countEvents(accountForCoEvent)
			for k in tmp.keys():
				if not devent.has_key(k):
					devent[k] = 0
				devent[k] += tmp[k]

		return devent