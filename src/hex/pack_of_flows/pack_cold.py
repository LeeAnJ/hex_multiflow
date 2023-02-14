from src.hex.pack_of_flows.pack_abc import HexPackOfFlows
from src.hex.single_flow.flow_cold import HexFlowCold


class HexPackOfColdFlows(HexPackOfFlows):
    def __init__(self, flows: list[HexFlowCold], dlt_eps: float):
        super().__init__(flows, dlt_eps)

    # ------------------------------------------------------------------------------------------------------------------
    # METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--------
    # ------------------------------------------------------------------------------------------------------------------
    def thsdq_profiling(self, i_operating_mode=0):
        super().thsdq_profiling(i_operating_mode)

    def calc_dq_w_distribution(self) -> [[], float]:
        # dlt_q_w - с нарастающим эффектом, т.е. dlt_q_w[n]-это все полученное тепло
        n = self.flow[0].n_nodes
        dq_w = []
        for i in range(n):
            dq_w_section = 0
            for flow in self.flow:
                dq_w_section += flow.dq_w[i]
            dq_w.append(dq_w_section)
        dq_w_sum = dq_w[-1]
        return dq_w, dq_w_sum

    def calc_ds_wk_distribution(self) -> [[], float]:
        n = self.flow[0].n_nodes
        ds_wk = []
        for i in range(n):
            ds_wk_section = 0
            for flow in self.flow:
                ds_wk_section += flow.ds_wk[i]
            ds_wk.append(ds_wk_section)
        ds_wk_sum = sum(ds_wk)
        return ds_wk, ds_wk_sum

    def dtmin_inter_flows_is_ok(self, dt_min: float = 0.0) -> bool:
        # находим минимальную разность температур между ближайшими парами потоков в пакете
        # ВАЖНО: списки потоков для пакетов были отсортированы по убыванию t_in (от большей t_in к меньшей)
        # еще при их создании в ..._flows_inp
        return super().dtmin_inter_flows_is_ok(dt_min)

    def _set_sections_in_flows(self, n_sections: int, n_nodes: int) -> None:
        super()._set_sections_in_flows(n_sections, n_nodes)

    def _combine_flows_outdat_in_sets(self):
        return super()._combine_flows_outdat_in_sets()

    def get_flows_t_in_lim(self):
        # for cold pack: t_in_lim_k = min(t_in_k)
        return min([self.flow[i].inlet.t_k for i in range(self.n_flows)])

    def calc_flows_out_thsds_with_eps_qact(self, q_act_w) -> []:  # [FlowOutData]:
        return super().calc_flows_out_thsds_with_eps_qact(q_act_w)

    def from_flows_remove_needless_eps(self, _mins: [float], _maxs: [float]) -> None:
        super().from_flows_remove_needless_eps(_mins, _maxs)

    def set_eps_per_flows(self, q_act_w: float, eps_hex: float) -> None:
        super().set_eps_per_flows(q_act_w, eps_hex)

    def _set_flows(self, flows) -> [int, []]:
        return super()._set_flows(flows)

    def _calc_flows_total_q_w(self) -> float:
        return super()._calc_flows_total_q_w()

    def calc_flows_eps_min_eps_max_with_q_act_w(self, q_act_w: float) -> [[], []]:
        return super().calc_flows_eps_min_eps_max_with_q_act_w(q_act_w)

    def set_flows_eps_from_eps_min_to_eps_max(self):
        return super().set_flows_eps_from_eps_min_to_eps_max()

    def combine_flows_epss_in_eps_sets(self) -> []:
        return super().combine_flows_epss_in_eps_sets()

    def get_eps_min_max_from_pack_eps_combined(self) -> [[], []]:
        return super().get_eps_min_max_from_pack_eps_combined()

    # def calc_flows_out_thsds_with_eps_qact(self, q_act_w) -> None: # [FlowOutData]:
    #     super().calc_flows_out_thsds_with_eps_qact(q_act_w)
