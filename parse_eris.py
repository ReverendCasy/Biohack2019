import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
import os

from keras.models import Model,Sequential
from keras.layers import *
from keras.utils import to_categorical
from keras.optimizers import SGD,Adam

import pandas as pd
from sklearn import datasets, linear_model
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import mean_squared_error, r2_score

datka = pd.read_csv('neur_net.txt', sep='\t')

concatdf = pd.concat([datka.toxin.apply(lambda x: pd.Series(list(x))),
           datka[' receptor'].apply(lambda x: pd.Series(list(x)))], axis=1)

concatdf.columns = [f'l_{i}' for i in range(concatdf.shape[1])]



from sklearn.preprocessing import LabelBinarizer
lb = LabelBinarizer()

X = pd.get_dummies(concatdf)
y = datka[' score']

from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4, random_state=101)

from sklearn.linear_model import LinearRegression
lm = LinearRegression()
lm.fit(X_train,y_train)

predictions = lm.predict(X_test)

plt.scatter(y_test,predictions);



