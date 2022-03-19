import Navbar from '../../components/Navbar/Navbar';
import Hero from '../../components/Hero/Hero';
import GetStarted from '../../components/GetStarted/GetStarted';
import styled from "styled-components";
import { Link } from 'react-router-dom';

function Home() {
    return (
        <backGround>
            <Navbar />
            <Hero />
            <GetStarted />
        </backGround>
    );
}

export default Home;

const backGround = styled.div`
    height: 100vh;
    background: #000;
`