from enum import Enum

class State(Enum):
	new="new"
	wait_for_registration="waiting for registration"
	wait_for_data="waiting for data"
	pipeline="pipeline"
	reprocess="reprocess"
	complete="complete"
	error="error"
