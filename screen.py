from kivy.app import App
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.screenmanager import ScreenManager, Screen
from kivy.core.window import Window
from kivy.config import Config
from kivy.factory import Factory
from kivy.properties import ObjectProperty
from kivy.uix.popup import Popup
from kivy.uix.checkbox import CheckBox
import webbrowser
from Bio.Seq import *
from Bio.Restriction import *
from Bio import SeqIO
import numpy
import json
import back as b

#Definir o tamanho que a tela deve ter e mandar ela abrir em full screen
Config.set('graphics', 'resizable', '0') #0 being off 1 being on as in true/false
Config.set('graphics', 'width', '1366')
Config.set('graphics', 'height', '700')
Config.write()
class Gerenciador(ScreenManager):
    pass

class Pre(Screen):
    tarefas = []
    path = ''

    def saveData(self,*args):
        with open(self.path+'data.json','w') as data:
            json.dump(self.tarefas,data)

    def addWidget(self):
        texto = self.ids.texto.text
        self.ids.texto.text = ''
        self.tarefas.append(texto)
        self.saveData()

class Aut(Screen):
    def on_pre_enter(self):
        Window.bind(on_keyboard=self.voltar)

    def voltar(self,window,key,*args):
        if key == 27:
            App.get_running_app().root.current = 'menu0'
            return True

    def on_pre_leave(self):
        Window.unbind(on_keyboard=self.voltar)
    def optimize(self):
        b.optim()
    def cod(self,cod):
        b.optimcod(cod)
    def on_checkbox_active(checkbox, value):
        if value:
            print('The checkbox', checkbox, 'is active')
        else:
            print('The checkbox', checkbox, 'is inactive')

    def tiraenzima(self,ecori,spei,psti,xbai):
    
        if ecori:
            b.ecoriswitch()
        if spei:
            b.speiswitch()
        if psti:
            b.pstiswitch()
        if xbai:
            b.xbaiswitch()

class Bases(Screen):
    def on_pre_enter(self):
        Window.bind(on_keyboard=self.voltar)

    def voltar(self,window,key,*args):
        if key == 27:
            App.get_running_app().root.current = 'aut'
            return True

    def on_pre_leave(self):
        Window.unbind(on_keyboard=self.voltar)

    def val1(self):
        webbrowser.open("https://rna.urmc.rochester.edu/RNAstructureWeb/Servers/Predict1/Predict1.html")
    def salvafasta(self):
        b.criarfasta('Sequencia otimizada por GammaJr','Sequencia')
    def gerarna(self):
        b.hairpincheckseq()

class Screen(App):
    def build(self):
        return Gerenciador()

Screen().run()
