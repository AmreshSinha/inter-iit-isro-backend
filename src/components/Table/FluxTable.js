import React from "react";
import styled from "styled-components";

function Table(props) {
  console.log(props);
  // console.log(Object.keys(props.props.Peak_Flux__x_).length);
  var indents = [];
  for (
    var i = 0;
    i < Object.keys(props.props.Peak_Flux__x_).length;
    i++
  ) {
    indents.push(
      <TableWrapper>
        <BurstStartX>{props.props.Peak_Flux__x_[`${i}`]}</BurstStartX>
        <BurstStartY>{props.props.Peak_Flux__y_[`${i}`]}</BurstStartY>
        <BurstPeakX>{props.props.background_count_Flux_vs_Time[`${i}`]}</BurstPeakX>
        <BurstPeakY>{props.props.Classification_by_Flux_Peak[`${i}`]}</BurstPeakY>
      </TableWrapper>
    );
  }
  return indents;
}

export default Table;

const TableWrapper = styled.div`
  background: rgba(237, 237, 237, 0.2);
  padding: 26px 36px;
  height: fit-content;
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  align-items: center;
  gap: 63px;
  border-radius: 12px;
  margin-bottom: 12px;
`;

const BurstStartX = styled.p`
  font-size: 1rem;
`;

const BurstStartY = styled.p`
  font-size: 1rem;
`;

const BurstPeakX = styled.p`
  font-size: 1rem;
`;

const BurstPeakY = styled.p`
  font-size: 1rem;
`;

const BurstEndX = styled.p`
  font-size: 1rem;
`;

const BurstEndY = styled.p`
  font-size: 1rem;
`;

const TotalBurstTime = styled.p`
  font-size: 1rem;
`;

const RiseTime = styled.p`
  font-size: 1rem;
`;

const DecayTime = styled.p`
  font-size: 1rem;
`;
