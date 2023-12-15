from typing import Optional

from pydantic import BaseModel


class SpiderSystem(BaseModel):
    source_name: str
    pos_ra: float 
    pos_dec: float 
    