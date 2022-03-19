import styled from "styled-components";

function Hero() {
    return(
        <HeroWrapper>
            <Welcome>Welcome to</Welcome>
            <HeroText>ISRO'S Web-Based<br />automatic identification of<br />Solar Bursts in X-Ray Light Curves</HeroText>
        </HeroWrapper>
    )
}

export default Hero;

const HeroWrapper = styled.div`
    margin-top: 58px;
    margin-left: 120px;
    display: flex;
    flex-direction: column;
`;

const Welcome = styled.h1`
    font-size: 64px;
    color: #F6C96F;
`;

const HeroText = styled.h1`
    font-size: 64px;
`