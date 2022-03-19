import styled from "styled-components";

function GetStarted() {
    return(
        <GetStartedWrapper>
            <a href='/input' style={GetStartedButton}>
                Get Started →
            </a>
        </GetStartedWrapper>
    )
}

export default GetStarted;

const GetStartedWrapper = styled.div`
    margin-top: 113px;
    margin-left: 120px;
`

const GetStartedButton = {
    color: '#F6C96F',
    fontSize: '32px',
    lineHeight: '24px',
    textDecoration: 'underline',
}