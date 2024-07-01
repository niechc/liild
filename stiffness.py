


class DoubleLinear:
    epsilon = 1e-6
    def __init__(self,
                 ke: float,
                 kp: float,
                 xy: float) -> None:
        self.ke = ke
        self.kp = kp
        self.xy = xy
        def upper(x: float) -> float:
            return self.xy * (self.ke - self.kp) + self.kp * x
        def lower(x: float) -> float:
            return - self.xy * (self.ke - self.kp) + self.kp * x
        self.upper = upper
        self.lower = lower
    
    def __call__(self, x: float, fs: float, dx: float) -> tuple:
        # upper yeild
        if abs(self.upper(x) - fs) < DoubleLinear.epsilon:
            if dx >= 0.0:
                return dx * self.kp + fs, self.kp
            else:
                return dx * self.ke + fs, self.ke
        # lower yeild
        elif abs(self.upper(x) - fs) < DoubleLinear.epsilon:
            if dx <= 0.0:
                return dx * self.kp + fs, self.kp
            else:
                return dx * self.ke + fs, self.ke
            
        # not yeild yet
        else:
            if self.upper(x + dx) < fs + dx * self.ke:
                return self.upper(x + dx), self.kp
            elif self.lower(x + dx) > fs + dx * self.ke:
                return self.lower(x + dx), self.kp
            else:
                return dx * self.ke + fs, self.ke
            
 
class Linear:
    def __init__(self,
                 k: float) -> None:
        self.k = k
        
    def __call__(self, x: float, fs: float, dx: float) -> float:
        return dx * self.k + fs, self.k


if __name__ == '__main__':
    pass