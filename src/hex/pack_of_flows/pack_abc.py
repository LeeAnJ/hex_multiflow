# import sys
import sys
from abc import ABC, abstractmethod
from itertools import product


class HexPackOfFlows(ABC):
    # eps-пакета: self.eps - это список, составленный из элементов списков eps-потоков: self.flow[i].eps
    # [ [self.flow[0][i].eps, self.flow[1][j].eps], ..., [self.flow[0][k].eps, self.flow[1][l].eps] ]
    # при условии равенства эл-тов суммы каждого подсписка 1: sum([self.flow[0][i].eps, self.flow[1][j].eps]) = 1
    # sum([self.flow[0][k].eps, self.flow[1][l].eps]) = 1

    # flow: список экземпляров потоков (горячих или холодных), кот. были иниц. еще в hot/cold_flows_inp
    def __init__(self, flow: [], dlt_eps: float):
        # потоки пакета были иниц. в hex_flows_init; здесь мы их просто переприсваиваем списку пакета self.flow[i]
        # ВАЖНО: списки потоков для пакетов были отсортированы по убыванию t_in (от большей t_in к меньшей)
        # еще при их создании в ..._flows_inp
        self.n_flows, self.flow = self._set_flows(flow)

        # предельная температура для каждого потока пакета на входе (определяем по входным тем-рам потоков пакета):
        # эти тем-ры будут использованы при расчете интервалов qt-диагр.: n=(t_in_lim_hot - t_in_lim_cold)/dt
        # hot pack: t_in_lim_k = max(t_in_k); cold pack:  = min(t_in_k)
        self.t_in_lim_k = self.get_flows_t_in_lim()

        # УДАЛИТЬ_________________________________________
        # предельная температура для каждого потока пакета на выходе (определяем по тем-рам встреч. потоков на входе):
        # hot pack: t_out_lim_k = min(t_in_hot_k); cold pack:  = max(t_in_cold_k)
        # self.t_out_lim_k = self.flow[0].outlet_lim.t_k
        # УДАЛИТЬ_________________________________________

        # sum_of_q_lim_w = sum(q_lim_w)
        self.sum_of_q_lim_w = self._calc_flows_total_q_w()

        # шаг для расчета векторов eps=q_w[i]/q_act_w, где i - индекс потока в пучке
        self.dlt_eps = dlt_eps

        # каждый поток имеет свой массив eps: self.flows[i].eps; их максимальные значения храним в массиве self.eps_max
        # очевидно, что eps_min = self.dlt_eps (поток не может иметь нулевую долю: eps = 0)
        self.eps_max = []
        self.eps_min = []

        self.eps_sets = []
        self.has_min_q_w = False  # this pack has minimal sum(q[i]), i E [1..n_flows]

        self.flows_outlet_dat = None
        self.flows_outlet_dat_in_sets = None  # see Hex func. calc_qt
        self.ds_wk_sect = []

        # изменение кол-ва тепла dq_w[i] = h[i]-h[0] и энтропии ds_wk[i]=s_wk[i]-s_wk[i-1] по сеч.ПАКЕТА
        # и суммарное изменение dq_w и ds_wk для всего пакета
        self.dq_w = []
        self.dq_w_sum = None
        self.ds_wk = []
        self.ds_wk_sum = None

    # ------------------------------------------------------------------------------------------------------------------
    # METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--------
    # ------------------------------------------------------------------------------------------------------------------
    @abstractmethod
    def thsdq_profiling(self, i_operating_mode=0):
        """теплообмен между потоками внутри одного пакета (горячего или холодного) не подразумевается;
           неважно пересекаются ли t(x) потоков в пределах одного пакета;
        """

        # для "пассивного" пакета возможен лишь ОДИН рабочий режим, т.е. i_operating_mode=0
        # для "активного" пакета рабочих режимов должно быть > 1.
        # т.к. у меня нет доч.класса для пассивного/активного пакета, а лишь для горячего/холодного,
        # то я вынужден оставить эту проверку активность-пассивность пакета
        if self.has_min_q_w and i_operating_mode != 0:
            print(f'for passive pack argument i_operating_mode != 0 in func. thsdq_profiling')
            sys.exit('critical error -> program terminated')

        # для заданного рабочего режима перебираем данные всех потоков пакета
        for i_flow in range(self.n_flows):
            # извлекаем параметры потоков на выходе
            h_out = self.flows_outlet_dat_in_sets[i_operating_mode][i_flow].h_jmol
            t_out = self.flows_outlet_dat_in_sets[i_operating_mode][i_flow].t_k
            s_out = self.flows_outlet_dat_in_sets[i_operating_mode][i_flow].s_jmolk
            # расчет распределений т/д параметров в потоках пакета
            self.flow[i_flow].t_k, \
                self.flow[i_flow].h_jmol, \
                self.flow[i_flow].s_jmolk, \
                self.flow[i_flow].ds_wk, \
                self.flow[i_flow].dq_w, \
                self.flow[i_flow].dq_w_sum = \
                self.flow[i_flow].calc_t_h_s_ds_distribution(t_out=t_out, h_out=h_out, s_out=s_out)
        # расчет распределения тепла и производство энтропии по секциям ПАКЕТА в целом (выше было по отдельным потокам)
        self.dq_w, self.dq_w_sum = self.calc_dq_w_distribution()
        self.ds_wk, self.ds_wk_sum = self.calc_ds_wk_distribution()

    @abstractmethod
    def calc_dq_w_distribution(self) -> [[], float]:
        pass

    @abstractmethod
    def calc_ds_wk_distribution(self) -> [[], float]:
        pass

    @abstractmethod
    def dtmin_inter_flows_is_ok(self, dt_min: float = 0.0) -> bool:
        # находим минимальную разность температур между ближайшими парами потоков в пакете
        # ВАЖНО: списки потоков для пакетов были отсортированы по убыванию t_in (от большей t_in к меньшей)
        # еще при их создании в ..._flows_inp
        is_ok = True
        if self.n_flows == 1:
            return is_ok
        for i in range(self.n_flows - 1):
            dlt_t = self.flow[i].t_k - self.flow[i + 1].t_k  # self.flow[i].t_k is np.array
            dt_min_inter_flows = min(dlt_t)
            if dt_min_inter_flows < dt_min:  # т.е. есть пара потоков, чьи температурные кривые пересеклись
                return not is_ok  # или сблизились больше допустимого
            return is_ok

    @abstractmethod
    def _set_sections_in_flows(self, n_sections: int, n_nodes: int) -> None:
        for flow in self.flow:
            flow.section_flow_for_qt_calc(n_sections=n_sections, n_nodes=n_nodes)

    @abstractmethod
    def _combine_flows_outdat_in_sets(self):

        data_sets = []

        for eps_set in self.eps_sets:
            # data_set_with_ds_sum = []  # наборы данных вместе с производством энтропии для последующ. сортировки
            data_set = []
            # ds_sum = 0.0
            for i_flow, eps in enumerate(eps_set):
                i_eps = self.flow[i_flow].eps.index(eps)

                data_set.append(self.flows_outlet_dat[i_flow][i_eps])

                # ds_sum += self.flows_outlet_dat[i_flow][i_eps].ds_wk

            # объединяем кортеж данных = f(eps) и суммарную энтропию в один список
            data_sets.append(tuple(data_set))
            # data_set_with_ds_sum.append(ds_sum)

            # добавим раб.режим: набор выход.данных как f(eps) и его суммар.энтропию я общий набор раб.режимов пакета
            # data_sets_with_ds_sum.append(data_set_with_ds_sum)
        return tuple(data_sets)

    @abstractmethod
    def get_flows_t_in_lim(self):
        # for hot pack: t_in_lim_k = max(t_in_k)
        # for cold pack: t_in_lim_k = min(t_in_k)
        pass

    @abstractmethod
    def calc_flows_out_thsds_with_eps_qact(self, q_act_w) -> ():
        # каждый поток пакета имеет свой кортеж eps: (eps_0, eps_1,...,eps_N ). для каждого eps_i проводим расчет
        # соответствующих выходных данных и сохраняем их в виде кортежа данных: (dat_0, dat_1,...,dat_N )
        # кортежи данных для каждого потока сохр. в кортеже: ( (dat_0, dat_1,...,dat_N), (dat_0, dat_1,...,dat_N ),...)
        # т.е. данные для i-го потока пакета и его j-го eps будем в дальнейшем находить по индексу dat[i][j]
        flow_outlet_dat_per_pack = []
        for flow in self.flow:
            outlet_dat_per_eps = []
            for eps in flow.eps:
                # thsdsdt_out: FlowOutData = (t, h, s, ds, dt) calc. at the flow's outlet port
                thsdsdt_out = flow.calc_out_thsdsdt_with_eps_qact(eps=eps, q_act_w=q_act_w)
                outlet_dat_per_eps.append(thsdsdt_out)  # формируем для каждого потока свой набор out dat = f(eps)
            # полученные списки для каждогого потока (dat, dat, dat) коррелируют со списком потока (eps, eps, eps)
            # списки выход.данных потоков отдного пакета собираем в один список (ну и конверт. в кортеж)
            flow_outlet_dat_per_pack.append(tuple(outlet_dat_per_eps))
        return tuple(flow_outlet_dat_per_pack)

    @abstractmethod
    def from_flows_remove_needless_eps(self, _mins: [float], _maxs: [float]) -> None:
        for i in range(self.n_flows):
            self.flow[i].remove_needless_eps(_min=_mins[i], _max=_maxs[i])

    @abstractmethod
    def set_eps_per_flows(self, q_act_w: float, eps_hex: float) -> None:

        if self.n_flows == 1:
            self.flow[0].eps = tuple([1.0])  # [1]
            return

        # self.n_flows > 1
        for i in range(self.n_flows):
            # "passive" pack отдает/получает все тепло, кот. он располагает -> выход.тем-ры его потоков стремятся
            # к своим предельным значениям (т.е. они "заданы"), а это значит, что теплосодерж. потоков (и их отношение)
            # тоже заданы однозначно. eps[i]=q[i]/q_act для каждого из потоков пассив.пакета задано однозначно!!!
            if self.has_min_q_w:

                self.flow[i].eps = self.flow[i].set_eps_for_flows_in_passive_pack(eps_hex=eps_hex, q_act_w=q_act_w)

                """ я предельные значения t_out_lim для потоков определяются в hex_flows_init. предельные значения
                h, s рассчитываются в __init__ flow: t=t_lim,h,q = f(t_lim, p).
                проблема в том, что в дальнейшем расчете t,h,q = f(q_act, eps) при eps_hex = 1 я опять "попаду"
                в "предельные" состояния, но ошибка расчетов, скажем, тем-ры в 0.13 К приводит к превышению
                ранее рассчитанных предельных значений и остановке программы.
                чтобы этого избежать, я хочу здесь пересчитать эти предельные значения: t,h,q = f(q_act,
                eps_max = eps[-1]), чтобы потом уже сравниваться с ними.
                это будет "избыточным" кодом, но я пока не знаю, как эту проблему решить иначе """

                if abs(eps_hex - 1.0) < 0.00001:  # eps_hex == 1.0
                    # оставляю неизменными: self.flow[i].outlet_lim.t_k
                    # оставляю неизменными: self.flow[i].outlet_lim.h_jmol
                    # оставляю неизменными: self.flow[i].outlet_lim.s_jmolk
                    self.flow[i].ds_lim_wk = (self.flow[i].outlet_lim.s_jmolk - self.flow[i].inlet.s_jmolk) * self.flow[
                        i].m_mols
                    self.flow[i].q_lim_jmol = self.flow[i].outlet_lim.h_jmol - self.flow[i].inlet.h_jmol
                    self.flow[i].q_lim_w = self.flow[i].q_lim_jmol * self.flow[i].m_mols
                else:  # eps_hex < 1.0
                    # print('abs(eps_hex-1.0) = ', abs(eps_hex-1.0))
                    # thsds_out = [FlowOutData] : t, p (None), h, s, ds
                    thsds_out, self.flow[i].q_lim_jmol, self.flow[i].q_lim_w = \
                        self.flow[i].calc_out_thsds_limiting_with_eps_qact(eps=self.flow[i].eps[0],
                                                                           eps_hex=eps_hex,
                                                                           q_act_w=q_act_w)
                    self.flow[i].outlet_lim.t_k = thsds_out.t_k  # temperature,  [K];
                    self.flow[i].outlet_lim.h_jmol = thsds_out.h_jmol
                    self.flow[i].outlet_lim.s_jmolk = thsds_out.s_jmolk
                    # не хочу "привязывать" ds_lim_wk к ... .outlet_lim - это св-во потока, а не просто выход. порта
                    self.flow[i].ds_lim_wk = thsds_out.ds_wk

            # "active" pack
            else:
                self.flow[i].eps = self.flow[i].set_eps_for_flows_in_active_pack(n_flows_in_pack=self.n_flows,
                                                                                 dlt_eps=self.dlt_eps,
                                                                                 q_act_w=q_act_w)

    @abstractmethod
    def _set_flows(self, flows_inp_dat) -> [int, []]:
        # set flows in the pack
        n_flows = len(flows_inp_dat)  # number of flows in the hot/cold pack
        flow = []
        for flow_inp_dat in flows_inp_dat:
            flow.append(flow_inp_dat)  # , t_out_lim_k=self.t_out_lim_k)
        return n_flows, flow

    @abstractmethod
    def _calc_flows_total_q_w(self) -> float:
        q_total_w = 0.0
        for i in range(self.n_flows):
            q_total_w += self.flow[i].q_lim_w  # self.flow_hot[i].q_max_w < 0  self.flow_cold[i].q_max_w > 0
        return q_total_w

    @abstractmethod
    def calc_flows_eps_min_eps_max_with_q_act_w(self, q_act_w: float) -> [[], []]:
        q_w_min = q_act_w * self.dlt_eps
        q_w_max = q_act_w - q_w_min * (self.n_flows - 1)
        # self.flow[i].q_w рассчитан ранее, исходя из предельной тем-ры потока на выходе
        eps_min = [self.dlt_eps] * self.n_flows
        eps_max = [round(abs(self.flow[i].q_lim_w / q_w_max), 4) for i in range(self.n_flows)]
        return eps_min, eps_max

    @abstractmethod
    def set_flows_eps_from_eps_min_to_eps_max(self) -> None:
        for i in range(self.n_flows):
            self.flow[i].set_eps(eps_min=self.eps_min[i], eps_max=self.eps_max[i], dlt_eps=self.dlt_eps)

    @abstractmethod
    # собрать все возможные комбинации eps потоков текущего пакета, которые дают в сумме 1.0, в список с эл-ми tuple,
    # где каждый элемент списка - tuple - это набор eps, дающий распределение q_act_w по потокам текущего пакета
    def combine_flows_epss_in_eps_sets(self) -> []:
        eps_in_one_list = []
        for i in range(self.n_flows):
            eps_in_one_list.append(self.flow[i].eps)

        eps_combined = []
        for item in product(*eps_in_one_list):
            if sum(item) == 1.0:
                eps_combined.append(item)

        return eps_combined

    @abstractmethod
    def get_eps_min_max_from_pack_eps_combined(self) -> [[], []]:
        # после того как ерs потоков одного пакета были "скомбинированы" (см. self.combine_flows_epss_in_eps_sets),
        # может оказаться, что меньшие значения eps потоков оказались невостребованными (если при расчете
        # self.pack_... .set_eps_per_flows(q_act_w=self.q_act_w, eps_hex=self.eps_hex) актуальные значения
        # макс. eps=f(q_act) оказались ниже изначально заданных eps=f(dlt_eps, n_flows)). это следует из
        # соблюдения условия sum(eps)=1.0. Тогда следует убрать эти "ненужные" значения eps из списков eps
        # для потоков, чтобы в дальнейшем не тратить время на расчеты вых.параметров потоков с этими eps
        eps_min = []
        eps_max = []
        for i in range(self.n_flows):
            eps_min.append(min(self.eps_sets, key=lambda tup: tup[i])[i])
            eps_max.append(max(self.eps_sets, key=lambda tup: tup[i])[i])
        return eps_min, eps_max
