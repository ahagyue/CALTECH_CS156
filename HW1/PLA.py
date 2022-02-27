import numpy as np
from abc import *

#학습 데이터 생성
data = np.array([[1, 0], [1, 1], [0, 1], [0, 0]])
label = np.array([1, 0, 1, 0])

#함수 추상 클래스
class function(metaclass = ABCMeta):
    @abstractmethod
    def forward(self):
        pass
    
    @abstractmethod
    def backward(self):
        pass

class relu(function):
    def __init__(self):
        self.neg = None
    
    def forward(self, x):
        self.neg = x<0
        x[self.neg] = 0
        return x
    
    def backward(self, out):
        out[self.neg] = 0
        return out

class softmax(function):
    def __init__(self):
        self.y = None
        
    def forward(self, x):
        self.y = 1/(1 + np.exp(-x))
        return self.y
    
    def backward(self, x):
        
