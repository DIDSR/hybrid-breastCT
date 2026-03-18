from abc import ABC, abstractmethod


class BaseSystem(ABC):
    def __init__(self, cfg: dict):
        self.cfg = cfg

    @abstractmethod
    def validate(self) -> None:
        pass

    @abstractmethod
    def summary(self) -> dict:
        pass
