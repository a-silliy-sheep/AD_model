
#region---module loading-----------------------------------------------
from sklearn.ensemble import RandomForestClassifier 
from sklearn.metrics import confusion_matrix, classification_report, accuracy_score
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import train_test_split
import numpy as np
import pandas as pd
from sklearn.tree import export_graphviz
import random, joblib
import streamlit as st

    # Set the dataframe output alignment:
pd.set_option('display.unicode.ambiguous_as_wide', True)
pd.set_option('display.unicode.east_asian_width', True)
pd.set_option('display.width', 180)
# endregion

#region---structure function-------------------------------------------------
def printDf(df):
    print(f"\033[1;34ma. 列名：\033[1;0m", f"{df.columns}", 
            f" \033[1;34m\nb. 各列数据类型：\033[1;0m", f"\n{df.dtypes}",
            f"\033[1;34m\nc. 维度: \033[1;0m", f"{df.shape}",
            f"\033[1;34m\nd. 前5行3列数据为: \033[1;0m", f"\n{df.iloc[0:5, 0:3]}")

def get_compositionDf(df):  # View composition ratios of categorical variables
    for i in df.columns:
        dat = df[i].value_counts()
        dat_ind = dat.index
        dat_ind.name = 'Factor'
        dat2 = pd.DataFrame({'count': dat.values, 'ratio (%)':np.round(100*dat.values/sum(dat.values), 1)}, index= dat_ind)
        print(f'\033[1;34m {i}: \n\033[1;0m', dat2, '\n')

# endregion ----------

#region ---------data loading --------------------------


df = pd.read_excel(r"D:\Code_Program\python\模型部署\atopic_dermatitis_250122\元数据\df_ml_python.xlsx")
X_train, X_test, y_train, y_test = train_test_split(df.iloc[:, 2:], df.iloc[:, 0], test_size=0.2, random_state=2025)



# endregion


# region----Model training and evaluation----------------------------------------
#           clf = RandomForestClassifier(n_estimators=100, random_state=2024)
#           clf.fit(X_train, y_train)
clf =  joblib.load(r".\RF.pkl")
y_pred = clf.predict(X_test)
y_pred2 = (clf.predict_proba(X_test)[:, 1] > 0.3).astype(int)
#       
print("\033[1;34mConfusion Matrix:\n\033[1;0m", pd.DataFrame(data=confusion_matrix(y_test, y_pred2),
            index=pd.Series(['0', '1'], name='Gold'),
            columns=pd.Series(['0', '1'], name='Diagnose')))  # 输出混淆矩阵
#     print("\033[1;34mClassification Report:\n\033[1;0m", classification_report(y_test, y_pred2, target_names=['0', '1']))  # 输出混淆矩阵衍生的各指标
#     print("Accuracy:\n", accuracy_score(y_test, y_pred))
#     print(clf.score(X_test, y_test))
#     clf.feature_names_in_
#     importances = clf.feature_importances_  # 计算特征重要性
#     print(importances)


# region-----streamlit Deployment of the online version of the forecasting tool---------------------

st.title('A simple tool to predict :blue[atopic dermatitis] among 2-8 year old preschool child',)

# number = st.sidebar.slider('选择一个数字', 0, 100, 50)  # 数据展示条
# data = pd.DataFrame({'a':np.random.randn(10), 'b':np.random.randn(10)})

# region----
#           if st.button('显示消息'):  # 按钮
#               st.write('Streamlit 是真的酷！')
#           if st.checkbox('显示图表'):  # 复选框
#               st.line_chart([0, 1, 2, 3, 4])
#           
#           genre = st.radio(
#             "你最喜欢哪种类型的音乐？",
#             ('流行', '摇滚', '爵士')
#           )
#           st.write(f'你的选择是：{genre}')
#           
#           age = st.slider('你的年龄：', 0, 130, 25)
#           st.write('我', age, '岁')
#           
#           # 初始化一个计数器
#           if 'count' not in st.session_state:
#               st.session_state.count = 0
#               
#           # 创建一个增加计数的按钮
#           if st.button('增加'):
#               st.session_state.count += 1
#           
#           st.write('计数器', st.session_state.count)
# endregion------

question_col_yesNo = ['Has father ever suffered from asthma?',
                'Has father ever suffered from allergic rhinitis?',
                'Has father ever suffered from atopic dermatitis?',
                'Has mother ever suffered from allergic rhinitis?',
                'Has mother ever suffered from atopic dermatitis?',
                'Renovating the dwelling before mother"s pregnancy?',
                'Renovating the dwelling during mother"s pregnancy?',
                'Is there any mold in the_dwelling during child first year of life?',
                'Does father smoking in the dwelling during child first year of life?',
                'Does the child have older brothers and sisters?',
                'Is the child birth weight less than 2500g?']
question_col_multi = ['Months of exclusive breastfeeding:',
                'Times of antibiotic therapy during child first year of life:']
st.sidebar.write('Blank in these questions:\n')

X_testStream = [1 if st.sidebar.selectbox(i, ('yes', 'no')) == 'yes' else 0 for i in question_col_yesNo]
last_q1 = st.sidebar.selectbox(question_col_multi[0], ('no', '1~3 months', '4 months or above')) 
last_q2 = st.sidebar.selectbox(question_col_multi[1], ('no', '1~2 times', '3 times or above'))


X_testStream = X_testStream + [(1 if ('no' in  last_q1) else 2 if ('1~3 months' in last_q1) else 3)
                ] + [(1 if ('no' in  last_q2) else 2 if ('1~2 times' in last_q2) else 3)]


X_testStreamDf = pd.DataFrame(data=np.array([X_testStream]),
                                columns=pd.Series(X_train.columns, name='columns'))


st.empty()
st.write('Data maping: \n', X_testStreamDf)
st.empty()

if st.sidebar.button('Predict'):
    y_predSt = clf.predict_proba(X_testStreamDf)[:,1][0]
    y_predSt2 = 0.15+((0.9-0.15)/(0.541-0.003))*(y_predSt-0.003) # 原区间[0.003, 0.541]，转换后区间[0.15, 9]
    st.markdown('''Based on the information you provided, the random forest model predicts the probability
        of your child developing atopic dermatitis between the ages of 2-8 as follows:''')
    st.markdown(f":red[{100*round(y_predSt2, 3)}%]")

#   if st.button('Predict'):
#       y_predSt = clf.predict_proba(X_testStreamDf)[:,1][0]
#       y_predSt2 = y_predSt * 0.5/0.24 if y_predSt <= 0.24 else 0.5 + ((1-0.5)/(1-0.25))*(y_predSt - 0.25)
#       st.markdown('''Based on the information you provided, the random forest model predicts the probability
#           of your child developing atopic dermatitis between the ages of 2-8 as follows:''')
#       st.markdown(f":red[{round(y_predSt2, 3)}]")
#   
# endregion

