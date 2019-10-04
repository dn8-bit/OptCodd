from kivy.app import App
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.screenmanager import ScreenManager, Screen
from kivy.core.window import Window
from kivy.config import Config
from kivy.factory import Factory
from kivy.properties import ObjectProperty
from kivy.uix.popup import Popup
from kivy.uix.checkbox import CheckBox
from kivy.uix.floatlayout import FloatLayout
import webbrowser
from Bio.Seq import *
from Bio.Restriction import *
from Bio import SeqIO
import numpy
import json
import back as b
import os
#Definir o tamanho que a tela deve ter e mandar ela abrir em full screen
Config.set('graphics', 'resizable', '0') #0 being off 1 being on as in true/false
Config.set('graphics', 'width', '1366')
Config.set('graphics', 'height', '700')
Config.write()
class Gerenciador(ScreenManager):
    pass

class Informacao(FloatLayout):
    pass
class Rna(FloatLayout):
    def salvafasta(self):
        return str(Seq(''.join(json.load(open('data.json')))).transcribe())[:40]
class SaveDialog(FloatLayout):
    save = ObjectProperty(None)
    text_input = ObjectProperty(None)
    cancel = ObjectProperty(None)



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
    def btn(self):
        show_popup()

class Bases(Screen):
    text_input = ObjectProperty(None)
    def on_pre_enter(self):
        Window.bind(on_keyboard=self.voltar)

    def voltar(self,window,key,*args):
        if key == 27:
            App.get_running_app().root.current = 'aut'
            return True
    def btn(self):
        show_popup()
    def salvin(self):
        content = SaveDialog(save=self.save, cancel=self.dismiss_popup)
        self._popup = Popup(title="Salvar arquivo", content=content,
                            size_hint=(0.9, 0.9))
        self._popup.open()

    def on_pre_leave(self):
        Window.unbind(on_keyboard=self.voltar)

    def val1(self):
        webbrowser.open("https://rna.urmc.rochester.edu/RNAstructureWeb/Servers/Predict1/Predict1.html")
    def salvafasta(self):
        return str(Seq(''.join(json.load(open('data.json')))).transcribe())[:40]
    def gerarna(self):
        show = Rna()
        popupWindow = Popup(title="40 primeiros pares de base", content=show, size_hint=(None,None),size=(500,300))
        popupWindow.open()
    def save(self, path, filename):
        with open(os.path.join(path, filename), 'w') as stream:
            stream.write(''.join(json.load(open('data.json'))))

        self.dismiss_popup()
    def dismiss_popup(self):
        self._popup.dismiss()


def show_popup():
    show = Informacao()
    popupWindow = Popup(title="INFORMAÇÃO", content=show, size_hint=(None,None),size=(250,150))
    popupWindow.open()


class Screen(App):
    def build(self):
        return Gerenciador()

Screen().run()
