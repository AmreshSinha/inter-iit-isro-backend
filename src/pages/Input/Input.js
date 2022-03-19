import styled, { keyframes } from "styled-components";
import Navbar from "../../components/Navbar/Navbar";
import { useState, useEffect, useRef } from "react";
import axios from "axios";
import CloudComputing from './Octicons-cloud-upload.svg';

function Input() {
    const [uploadState, changeUploadState] = useState(false);
    let [selectedFile, setSelectedFile] = useState(null);
    function onFileChange(event) {
        setSelectedFile( selectedFile= event.target.files[0] );
    };
    function onFileUpload() {
        changeUploadState(!uploadState);
        // Create an object of formData
        const formData = new FormData();

        // Update the formData object
        formData.append(
            'imgfile',
            selectedFile
        )

        // Details of the uploaded file
        console.log(selectedFile);

        // Sending data to backend
        // axios.post('http://172.20.10.2:8080/api/upload', formData).then(res => {
        //     changeUploadState(uploadState=true);
        //     console.log(uploadState)
        // })
        fetch('http://localhost:8080/api/upload', 
            {
                method: 'POST',
                body: formData
            }
        )
        .then(res => res.json())
        .then(res => {if (res.status == 'ok') {
            window.location = '/prediction'
        }})
        .catch(err => console.error('Error:', err))
    }
    function fileData() {
        if (selectedFile) {
            return (
                <h1>Uploaded!<br />Details: {selectedFile}</h1>
            )
        } else {
            return (
                console.log('No file selected')
            )
        }
    }
    return (
        <>
            <Navbar />
            <InputWrapper>
                <InputData>Input Data</InputData>
                <InputBoxWrapper>
                    <InputBox type='file' onChange={onFileChange} />
                    <CloudIcon src={CloudComputing}/>
                </InputBoxWrapper>
                {!uploadState ? <IconButtonWrapper>
                    <InputButton onClick={onFileUpload}>Upload</InputButton>
                </IconButtonWrapper> : <IconButtonWrapper><InputButton disabled style={{cursor: "none"}}>Uploading...</InputButton></IconButtonWrapper>}
            </InputWrapper>
            <InputInfoWrapper>
                File Format Support: FITS, XLS, ASCII
            </InputInfoWrapper>
        </>
    )
}

export default Input;

const InputWrapper = styled.div`
    width: fit-content;
    margin-left: 120px;
    background: rgba(237, 237, 237, 0.2);
    height: fit-content;
    display:  flex;
    flex-direction: column;
    gap: 65px;
    padding: 30px;
    border-radius: 12px;
    margin-bottom: 35px;
`
const InputData = styled.div`
    font-size: 32px;
    color: #F6C96F;
`

const InputBoxWrapper = styled.div`
    background: rgba(237, 237, 237, 0.2);
    border-radius: 8px;
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    align-items:  center;
    width: 469px;
    height: 64px;
    padding-left: 16px;
    padding-right: 16px;
`;

const CloudIcon = styled.img`
    width: 14px;
    height: 17px;
`

const InputBox = styled.input`
    width: 469px;
    border-radius: 8px;
    color: #fff;
    font-size: 16px;
`

const IconButtonWrapper = styled.div`
    width: 100%;
    display: flex;
    flex-direction: row;
    justify-content: flex-end;
`

const InputButton = styled.button`
    font-family: inherit;
    background: #F6C96F;
    border-radius: 8px;
    border: none;
    width: 86px;
    height: 56px;
    cursor: pointer;
    &:hover {
        border: 2px solid #fff;
    }
`

const InputInfoWrapper = styled.div`
    width: 501px;
    background: rgba(237, 237, 237, 0.2);
    margin-left: 120px;
    padding: 30px;
    border-radius: 12px;
`;

const OverlaySpinner = styled.div`
    margin: 0;
    padding: 0;
    position: fixed;
    z-index: 5;
    width: 100vw;
    height: 100vh;
    background: rgba(237, 237, 237, 0.2);
`

const spin = keyframes`
    to { 
        transform: rotate(360deg);
    }
`

const Spinner = styled.div`
    display: inline-block;
    width: 50px;
    height: 50px;
    border: 3px solid rgba(255,255,255,.3);
    border-radius: 50%;
    border-top-color: #fff;
    animation: ${spin} 1s ease-in-out infinite;
    /* -webkit-animation: spin 1s ease-in-out infinite; */
`